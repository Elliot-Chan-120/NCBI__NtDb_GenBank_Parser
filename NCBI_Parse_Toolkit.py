from Bio import Entrez, SeqIO
from pathlib import Path
import plotly.graph_objects as objects


class NCBIParser:
    """
    :param: accession_no
    :return: parsed genbank record -> extract key info like identifiers, seq details and feature annotations
    """
    def __init__(self, accession_no, email):
        self.accession_no = accession_no
        Entrez.email = email
        try:
            with Entrez.efetch(
                db="nucleotide",
                id=self.accession_no,
                rettype="gb",  # request genbank format
                retmode="text",  # request text mode
                seq_start=1  # ensure full handle
            ) as self.handle:
                self.record = SeqIO.read(self.handle, "genbank")
        except Exception as e:
            print(f"Failed to fetch record {str(e)}")

        # class lists to contain parsed data
        # can reuse them in (future) functions
        self.feature_data = []
        self.sequence_data = []

    def log_results(self):
        """Will add this in later"""
        pass

    def parse_seq_info(self):
        """Returns list of basic sequence information"""
        try:
            sequence_info = {
                #Identifiers
                'Accession': self.record.id,
                'Locus_name': self.record.name,
                'Version': self.record.annotations.get('sequence_version'),

                #Sequence Details
                'Sequence_length': len(self.record.seq),
                'Description': self.record.description,
                'Organism': self.record.annotations.get('organism'),
                'Taxonomy': self.record.annotations.get('taxonomy', []),
                'Full_sequence': str(self.record.seq),
            }
            self.sequence_data.append(sequence_info)
            return sequence_info
        except Exception as e:
            print(f"Failed to obtain sequence info {str(e)}")

    def parse_feature_info(self):
        """Returns list of dictionaries of all gene features"""
        self.feature_data = []

        try:
            for feature in self.record.features:
                if feature.type == "source":
                    source_data = {
                        'Type': feature.type.upper(),
                        'Location': feature.location,
                        'Start': int(feature.location.start),
                        'End': int(feature.location.end),
                        'Organism info': f"{feature.qualifiers.get('organism', [''])[0]}, "
                                         f"plasmid {feature.qualifiers.get('plasmid', [''][0])}",
                        'Molecular Type': feature.qualifiers.get('mol_type', [''])[0],
                    }
                    self.feature_data.append(source_data)
                elif feature.type == "gene":
                    gene_data = {
                        'Type': feature.type.upper(),
                        'Start': int(feature.location.start),
                        'End': int(feature.location.end),
                        'Gene': feature.qualifiers.get('gene', ['N/A'])[0],
                        'Locus Tag': feature.qualifiers.get('locus_tag', ['N/A'])[0],
                        }
                    self.feature_data.append(gene_data)
                elif feature.type == "CDS":
                    cds_data = {
                        'Type': feature.type.upper(),
                        'Start': int(feature.location.start),
                        'End': int(feature.location.end),
                        'Gene': feature.qualifiers.get('gene', ['N/A'])[0],
                        'Locus Tag': feature.qualifiers.get('locus_tag', ['N/A'])[0],
                        'Codon Start': feature.qualifiers.get('codon_start', ['1'])[0],
                        'Protein ID': feature.qualifiers.get('protein_id', ['N/A'])[0],
                        'Protein Prod. Description': feature.qualifiers.get('product', ['N/A'])[0],
                        'Translated AA Chain': feature.qualifiers.get('translation', ['N/A'])[0],
                        'Notes': feature.qualifiers.get('note', ['N/A'])[0],
                    }
                    self.feature_data.append(cds_data)
                elif feature.type == "misc_feature":
                    misc_data = {
                        'Type': feature.type.upper(),
                        'Location': str(feature.location),
                        'Start': int(feature.location.start),
                        'End': int(feature.location.end),
                        'note': feature.qualifiers.get('note', ['N/A']),
                        'db_xref': feature.qualifiers.get('db_xref', ['N/A'])
                    }
                    self.feature_data.append(misc_data)

            return self.feature_data
        except Exception as e:
            print(f"Failed to obtain feature info {str(e)}")

    def output_seq_info(self):
        """Outputs sequence information"""
        if not self.sequence_data:
            self.parse_seq_info() # make sure the seq info is parsed before printing anything (don't want empty list)

        for key, value in self.sequence_data[0].items(): # access the first (and only dictionary in the list)
            print(f"{key}: {value}") # print out its values

    def output_feature_info(self):
        """Outputs all feature dictionaries: source, gene, cds"""
        if not self.feature_data:
            self.parse_feature_info()

        for x in range(len(self.feature_data)):
            print("=" * 50)
            print(f"Feature Dictionary {x + 1}")
            for feature_x, key in self.feature_data[x].items():
                print(f"{feature_x}: {key}")

    def batch_parse(self, accession_list):
        """Will add this in later: -> Parse Multiple Sequences at once"""
        pass

    def calculate_seq_stats(self):
        """
        Will add this in later: Calculate GC content, codon usage -  should probably be able to import from a separate class/project
        """
        pass

    def linear_gene_map(self):
        """
        Create linear visual representation of gene locations
        !! When coming across GenBank ID's with many feature dicts (40+), it can look clunky
        \ntry zooming in, trying to find a way to work around this
        """
        fig = objects.Figure()

        # draw main sequence line
        sequence_length = len(self.record.seq)

        fig.add_trace(objects.Scatter(
            x=[0, sequence_length],
            y=[0, 0],
            mode='lines',
            line=dict(color='black', width=2),
            name='Sequence',
            hoverinfo='skip'
        ))

        # Colours!
        color_map = {
            'cds': '#2ecc71',
            'gene': '#3498db',
            'source': '#95a5a6'
        }

        # Need varying heights for different dicts since they're probably going to overlap
        y_positions = {
            'SOURCE': 0.02,
            'GENE': -0.05,
            'CDS': -0.15,
        }

        for feature_dict in self.feature_data:
            f_type = feature_dict['Type'].upper()
            start = feature_dict['Start']
            end = feature_dict['End']

            # initialize y_pos
            y_pos = y_positions.get(f_type, 0)

            # Create hover texxt based on feature type
            hover_text = ""
            if f_type == "SOURCE":
                hover_text += (f"org_info {feature_dict['Organism info']} "
                               f"<br> mol_type: {feature_dict['Molecular Type']}")
            elif f_type == "GENE":
                hover_text += f"locus: {feature_dict['Locus Tag']}"
            elif f_type == "CDS":
                hover_text += (f"gene: {feature_dict['Gene']} "
                               f"<br>locus: {feature_dict['Locus Tag']} "
                               f"<br>protein_ID: {feature_dict['Protein ID']}"
                               f"<br>protein_desc: {feature_dict['Protein Prod. Description']}"
                               f"<br>First 20 AAs: {feature_dict['Translated AA Chain'][:20]}"
                               f"<br>notes: {feature_dict['Notes']}")

            # feature objects
            fig.add_trace(objects.Scatter(
                x=[start, start, end, end],             # Four x coords
                y=[0, -y_pos, -y_pos, 0],               # Four y coords for rect.
                fillcolor=color_map.get(f_type),
                line=dict(width=1),                     # No border
                name=f_type,                            # Name in legend
                text=hover_text,                        # Hover text created
                hoverinfo='text',                       # Show custom hover text
                showlegend=True                         # Show legend
            ))

            fig.add_annotation(
                x=(start+end) / 2,
                y=y_pos + 0.23,
                text=f"{f_type} <br> {start}-{end}",
                showarrow=False,
                font=dict(size=10)
            )

        # Update Layout
        fig.update_layout(
            title=f"Linear Gene Map: {self.accession_no}, {self.sequence_data[0]['Organism']}",
            xaxis=dict(
                title="Sequence Position (bp)",
            range=[-sequence_length * 0.05, sequence_length * 1.05  ]
            ),
            yaxis=dict(
                range=[-1, 1],
                showticklabels=False,
                showgrid=False,
                zeroline=False
            ),
            showlegend=True,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99
            ),
            hovermode='closest',
            plot_bgcolor='white'
        )

        fig.show()

    def NCBI_Parsed_File(self):
        """
        :return: file named "NCBI_parse_results" with parsed information
        """
        path = Path('NCBI_parse_results')
        content = ""
        # write sequence information
        content += "Sequence Information:\n"
        for k, v in self.parse_seq_info().items():
            content += f"{k}: {v} \n"

        # write feature information
        content += "\nFeature Information:\n"
        for number in range(len(self.parse_feature_info())):
            content += f"{"=" * 50}\nFeature Dictionary {number + 1}\n"
            for features, data in self.parse_feature_info()[number].items():
                content += f"{features}: {data}\n"

        # write out final file
        path.write_text(content)
