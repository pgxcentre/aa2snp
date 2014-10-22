#!/usr/bin/env python3

import logging
import argparse
import urllib.request
# urllib.request.urlopen(url[, data][, timeout])
import json
import re

# We will add build information after
ensembl_url = "rest.ensembl.org"


def ensembl_protein_to_genomic(protein, position):
    """Uses the Ensembl REST API to convert a protein position to the cDNA coordinate.

    :param protein: The Ensembl protein ID or symbol (e.g. ENSP00000288602, BRCA2).
    :param region: The amino acid position (e.g. 234).

    """

    if not protein.startswith("ENSP"):
        protein = symbol_lookup(protein)

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


def main():
    global ensembl_url

    args = parse_args()
    build = args.build
    if build == "GRCh38":
        build = "" # We assume Ensembl uses GRCh38 by default.
    else:
        ensembl_url = ".".join((build, ensembl_url))

    position = re.match("^[A-Z]([0-9]+)[A-Z]$", args.var)
    if position is None or not position.group(1):
        raise Exception("Invalid format for mutation notation: {}".format(
            args.var
        ))
    position = position.group(1)

    # Get potential genomic mappings.
    mappings = ensembl_protein_to_genomic(args.gene, position)

    # Query these mappings for SNPs.
    print("\t".join((
        "id",
        "chromosome",
        "position",
        "strand",
        "alleles"
    )))

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

        variants = variants_in_region(region)
        for var in variants:
            assert var["start"] == var["end"]
            print("\t".join([str(i) for i in (
                var["id"],
                var["seq_region_name"],
                var["start"],
                var["strand"],
                "/".join(var["alt_alleles"]))]
            ))


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
        help="The genomic build (default: GRCh37.",
        type=str,
        default="GRCh37",
        choices=("GRCh37", "GRCh38")
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()


