#!/usr/bin/env python3

import logging
import argparse
import urllib.request
import json
import re
import sys
from collections import defaultdict

# We will add build information after
ensembl_url = "rest.ensembl.org"


def ensembl_protein_to_genomic(protein, position):
    """Uses the Ensembl REST API to convert a protein position to the cDNA coordinate.

    :param protein: The Ensembl protein ID or symbol (e.g. ENSP00000288602, BRCA2).
    :param region: The amino acid position (e.g. 234).

    """

    if not protein.startswith("ENSP"):
        protein = symbol_lookup(protein)
        if protein is None:
            sys.exit("You can fix this problem by providing an Ensembl protein ID "
                "instead of a gene symbol.")

    region = "..".join((position, position))

    url = ("http://{}/map/translation/{id}/{region}?"
           "content-type=application/json&"
           "species=homo_sapiens")

    url = url.format(
        ensembl_url,
        id=protein,
        region=region
    )
    logging.debug("Queried url: " + url)

    with urllib.request.urlopen(url) as stream:
        res = json.loads(stream.read().decode("utf-8"))

    mappings = res.get("mappings")
    if mappings is None:
        logging.critical("Could not find mapping for mutation {} in protein "
            "{}".format(region, protein))
        return None
    
    return mappings


def symbol_lookup(symbol):
    """Converts a protein symbol to an Ensembl Protein ID (ENSP).

    :param symbol: A protein symbol (e.g. BRCA2).

    """

    url = ("http://{}/lookup/symbol/homo_sapiens/{}"
           "?content-type=application/json&expand=1")

    url = url.format(
        ensembl_url,
        symbol
    )
    logging.debug("Queried url: " + url)

    with urllib.request.urlopen(url) as stream:
        res = json.loads(stream.read().decode("utf-8"))

    # Look over transcripts.
    prot_id = None
    for transcript in res["Transcript"]:
        cur_prot_id = transcript["Translation"]["id"]
        if prot_id is None:
            prot_id = cur_prot_id
        elif prot_id != cur_prot_id:
            logging.critical("Ambiguous protein id for symbol {}".format(
                symbol
            ))
            return None

    return prot_id


def variants_in_region(region):
    """Retrieves the SNPs in a given region.

    :param region: Genomic region of the form C:123-456:S where C is the 
                   chromosome and S is the strand.

    """

    url = ("http://{}/overlap/region/human/{}?feature=variation"
           "&species=homo_sapiens"
           "&content-type=application/json")

    url = url.format(ensembl_url, region)
    logging.debug("Queried url: " + url)

    with urllib.request.urlopen(url) as stream:
        res = json.loads(stream.read().decode("utf-8"))

    return res


def variant_amino_changes(variant):
    """Retrieves variant consequences.

    :param variant: A JSON describing the variation.

    """

    # The region
    region = "{chrom}:{start}-{end}:{strand}".format(
        chrom=variant["seq_region_name"],
        start=variant["start"],
        end=variant["end"],
        strand=variant["strand"])

    url = ("http://{}/vep/homo_sapiens/region/{}/{}?"
           "content-type=application/json")

    # The amino changes
    amino_changes = set()

    # Cycling through the alternative alleles (the first allele is reference)
    for allele in variant["alt_alleles"][1:]:
        url = url.format(ensembl_url, region, allele)
        logging.debug("Queried url: " + url)

        with urllib.request.urlopen(url) as stream:
            results = json.loads(stream.read().decode("utf-8"))

            for res in results:
                for cons in res["transcript_consequences"]:
                    if "amino_acids" in cons:
                        amino_changes.add(cons["amino_acids"])

    return amino_changes


def main():
    global ensembl_url

    args = parse_args()
    build = args.build
    if build == "GRCh38":
        build = "" # We assume Ensembl uses GRCh38 by default.
    else:
        ensembl_url = ".".join((build, ensembl_url))

    matched = re.match("^([A-Z])([0-9]+)([A-Z])$", args.var)
    if matched is None:
        raise Exception("Invalid format for mutation notation: {}".format(
            args.var
        ))

    # Getting matched position and aminos
    position = matched.group(2)
    amino_ref = matched.group(1)
    amino_alt = matched.group(3)

    # Checking we have position and aminos
    if not position or not amino_ref or not amino_alt:
        raise Exception("Invalid format for mutation notation: {}".format(
            args.var
        ))

    # Get potential genomic mappings.
    mappings = ensembl_protein_to_genomic(args.gene, position)

    # Query these mappings for SNPs.
    print("\t".join((
        "#id",
        "chromosome",
        "position",
        "strand",
        "alleles",
        "amino_changes"
    )))

    # Amino acid problem
    has_amino_problem = False
    amino_problem = defaultdict(set)

    for mapping in mappings:
        # Build a region string.
        assert mapping["coord_system"] == "chromosome"
        assert build == mapping["assembly_name"]

        chrom = mapping["seq_region_name"]
        start = mapping["start"]
        end = mapping["end"]
        strand = mapping["strand"]

        region = "{}:{}-{}:{}".format(
            chrom, start, end, strand
        )
        # User can request that we print the genomic region corresponding
        # to the amino acid.
        if args.print_region:
            print("#Genomic mapping: {}".format(region))

        variants = variants_in_region(region)
        for var in variants:
            # Fetching the variant's consequences
            amino_changes = variant_amino_changes(var)

            for amino_change in amino_changes:
                var_amino_ref, var_amino_alt = amino_change.split("/")

                if var_amino_ref != amino_ref and var_amino_alt != amino_alt:
                    has_amino_problem = True
                    amino_problem[var["id"]].add(amino_change)

            # Printing the results
            print(
                var["id"],
                var["seq_region_name"],
                var["start"],
                var["strand"],
                "/".join(var["alt_alleles"]),
                ",".join(amino_changes),
                sep="\t"
            )

    # Are there any amino problem?
    if has_amino_problem:
        for marker, problems in amino_problem.items():
            logging.warning(
                "variant {} doesn't have right amino change {} "
                "instead of {}".format(
                    marker,
                    ",".join(problems),
                    amino_ref + "/" + amino_alt
                )
            )


def parse_args():
    """Parses command line arguments.

    --gene: The gene symbol.
    --var: The amino acid change AxxxB (A is WT and B is mutation).
    --build: Either GRCh37 or GRCh38 for now.

    """
    parser = argparse.ArgumentParser(
        description="Converts the amino acid change into a cDNA change."
    )

    parser.add_argument("-g", "--gene",
        help="The gene symbol for this mutation.",
        type=str,
        required=True,
    )

    parser.add_argument("-v", "--var",
        help=("The amino acid change (format is AxxxB where A is WT and B is " 
              "the new amino acid."),
        type=str,
        required=True,
    )
    
    parser.add_argument("-b", "--build",
        help="The genomic build (default: %(default)s).",
        type=str,
        default="GRCh37",
        choices=("GRCh37", "GRCh38")
    )

    parser.add_argument("-pr", "--print_region",
        help=("Print the genomic region corresponding to the given amino acid "
              "change."),
        action="store_true",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()


