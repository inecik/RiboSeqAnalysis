"""
Some functions previously written.
"""

import os
import numpy as np
import subprocess
import pandas as pd
from infrastructure.main import progress_bar


def gtf_cds_parser(gtf_path, verbose=False):
    """
    This is to convert gtf file into list of dictionary elements.
    :param verbose: Boolean. Report or not report the progress.
    :param gtf_path: String. Path of the gtf or gtf.gz file
    :return: List of dictionaries for each entry
    """
    # Column names for GTF file. Source: https://www.ensembl.org/info/website/upload/gff.html
    columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame"]
    assert os.path.splitext(gtf_path)[1] == ".gtf"
    with open(gtf_path, "r") as gtf_handle:  # Open and read the file
        gtf = gtf_handle.readlines()
    # Remove the titles, strip the trailing white spaces, split the line into cells
    gtf = [i.strip().split('\t') for i in gtf if not i.startswith('#')]
    output = list()
    filter_columns = ["protein_id", "transcript_id", "gene_name", "start", "end",
                      "strand", "frame", "seqname", "exon_number"]
    for ind, the_line in enumerate(gtf):  # Parse the entries
        # Parse the attributes column, column 9, of each line
        # noinspection PyTypeChecker
        entry = dict(  # Get the result as dictionary
            [[k, l.replace("\"", "")] for k, l in  # 4) Remove double quote from value
             [j.strip().split(" ", 1) for j in  # 3) For each attribute split key from value
              the_line[8].split(';')  # 1) Split first with ';'
              if j]])  # 2) See if element contains anything (last character is ';', to remove the artifact of it)
        entry.update(dict(zip(columns, the_line[:8])))  # Append dictionary with remaining information (Info in columns)
        # noinspection PyTypeChecker
        if entry["feature"] == "CDS":
            output.append([entry[ft] if ft in entry else np.nan for ft in filter_columns])
        if verbose and (ind % 100 == 0 or ind == len(gtf) - 1):  # Show the progress if 'verbose' is True
            progress_bar(ind, len(gtf) - 1)
    return pd.DataFrame(output, columns=filter_columns)


def download_gtf_refseq(temp_repo_dir, data_url=None):
    """
    Function is to download or find the file in temp_repo_dir.
    :param temp_repo_dir: Directory to download or find the file
    :param data_url: URL of the data.
    :return: String. Path of gtf or gtf.gz
    """
    if not data_url:
        data_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/" \
                   "GCF_000001405.38_GRCh38.p12_genomic.gtf.gz"
    gtf_file = os.path.basename(data_url)  # File name in the FTP server
    # Paths to download the database or look for the database if it already exists
    gtf_gz_path = os.path.join(temp_repo_dir, gtf_file)  # Compressed .gz format
    gtf_path = os.path.splitext(gtf_gz_path)[0]  # Not compressed
    if os.access(gtf_path, os.R_OK) or os.path.isfile(gtf_path):  # Check if the file exist
        return gtf_path  # Return the path if it exist
    else:  # Download otherwise
        subprocess.run((f"cd {temp_repo_dir}; "
                        f"curl -L -R -O {data_url}; "
                        f"gzip -d {gtf_file}"), shell=True)
        assert os.path.isfile(gtf_path) and os.access(gtf_path, os.R_OK)
        return gtf_path  # Return the compressed file
