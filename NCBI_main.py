from NCBI_Parse_Toolkit import *

result = NCBIParser("NC_005098.1", "examplename123@gmail.com")

result.parse_seq_info()

result.parse_feature_info()


# Output Functions - Displays information in terminal
result.output_seq_info()

result.output_feature_info()

result.NCBI_Parsed_File()

result.linear_gene_map()