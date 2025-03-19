import pandas as pd
from NCBI_Parse_Toolkit import NCBIParser


# Going to add gene cluster batch parsing
# Maybe a summary gene map comparing everything in one linear gene map? -> Might be too slow / unusable for big datasets

class BatchParse:
    """Given a list of accession IDs, it will batch parse them all and output their sequence information into a csv file"""
    def __init__(self, accession_list, email):
        self.accession_list = accession_list
        self.email = email

        # class lists containing batch parsed data to be reused in functions
        # easier than having to recall them from the previous class
        self.seqinfo_record = []
        self.featureinfo_record = []

    def all_info(self):
        """Fetch all sequence information"""
        print(self.accession_list)
        for accession in self.accession_list:
            acc_idx = NCBIParser(accession, self.email)
            acc_idx.parse_seq_info()
            self.seqinfo_record.append(acc_idx.sequence_data[0])


    def test_info(self):
        for x in range(len(self.accession_list)):
            for key, value in self.seqinfo_record[x].items():
                if key != 'Full_sequence':
                    print(f"{key}: {value}")


    def filecompare_seqinfo(self):
        """Output all sequence information into csv file for comparison"""
        process_data = []
        for idx in range(len(self.seqinfo_record)):
            add_info = {}
            for key, value in self.seqinfo_record[idx].items():
                if key != 'Full_sequence':
                    add_info[key] = value
            process_data.append(add_info)

        df = pd.DataFrame(process_data)
        df.to_csv('test.csv')

