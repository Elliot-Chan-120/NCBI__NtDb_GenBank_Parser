from NCBI_Parse_Toolkit import *
from Batch_Parser import BatchParse

result = NCBIParser("NC_001802", "example123@gmail.com")
result.parse_seq_info()
result.parse_feature_info()


def demo():
    result.output_seq_info()

    result.output_feature_info()

    result.NCBI_Parsed_File()

    result.linear_gene_map()

    result.gene_cluster(500)


accession_list = [
    "NC_000852",  # Example viral genome
    "NC_001802",  # HIV-1 genome
    "NC_003888",  # SARS-CoV genome
    "NC_002516",  # Pseudomonas aeruginosa genome
    "NC_000913"   # E. coli K-12 genome
]

accession = "NC_005098.1"

batch_result = BatchParse(accession_list, "elliotchan120@gmail.com")
batch_result.all_info()
batch_result.test_info()
batch_result.filecompare_seqinfo()


