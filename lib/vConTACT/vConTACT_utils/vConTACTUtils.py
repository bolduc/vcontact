import time
import os
import subprocess
import csv
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from string import Template
import pandas as pd

from KBaseReport.KBaseReportClient import KBaseReport


def log(message, prefix_newline=False):
    """
    Logging function, provides a hook to suppress or redirect log messages.
    """
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


html_template = Template("""<!DOCTYPE html>
<html lang="en">
  <head>
    <link href="https://netdna.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdn.datatables.net/1.10.19/css/jquery.dataTables.min.css" rel="stylesheet">
    <script src="https://code.jquery.com/jquery-3.3.1.js" type="text/javascript"></script>
    <script src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js" type="text/javascript"></script>
  </head>
  <body>
    <div class="container">
      <div>
        ${html_table}
      </div>
    </div>

    <script type="text/javascript">
      $$(document).ready(function() {
          $$('#my_id').DataTable();
      } );
      </script>
  </body>
</html>""")


class vConTACTUtils:

    def __init__(self, config):
        self.scratch = os.path.abspath(config['scratch'])

    def vcontact_help(self):
        command = "vcontact --help"
        self._run_command(command)

    def run_vcontact(self, params):

        #
        mappings = {
            'gene2genome': '--proteins-fp',
            'sequences': '--raw-proteins',
            'db': '--db',
            'pcs_mode': '--pcs-mode',
            'vcs_mode': '--vcs-mode',
            'blast_evalue': '--blast-evalue',
            'pc_max_overlap': '--max-overlap',
            'pc_penalty': '--penalty',
            'pc_haircut': '--haircut',
            'pc_inflation': '--pc-inflation',
            'vc_inflation': '--vc-inflation',
            'vc_density': '--min-density',
            'vc_min_size': '--min-size',
            'vc_max_overlap': '--vc-overlap',
            'vc_penalty': '--vc-penalty',
            'vc_haircut': '--vc-haircut',
            'merge_method': '--merge-method',
            'similarity': '--similarity',
            'seed_method': '--seed-method',
            'min_significance': '--sig',
            'max_significance': '--max-sig',
            'module_inflation': '--mod-inflation',
            'mod_significance': '--mod-sig',
            'module_min_shared': '--mod-shared-min',
            'link_significance': '--link-sig',
            'link_proportion': '--link-prop'
        }

        bool_args = ['optimize', 'permissive']

        # Should create build_command?
        command = 'vcontact --output-dir outdir'
        # Binaries
        command += ' --diamond-bin /usr/local/bin/diamond --c1-bin /usr/local/bin/cluster_one-1.0.jar'

        for param, cmd in mappings.items():
            command += ' {} {}'.format(cmd, params[param])

        self._run_command(command)

    def _run_command(self, command):
        """
        _run_command: run command and print result
        """

        log('Start executing command:\n{}'.format(command))
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output = pipe.communicate()[0]
        exitCode = pipe.returncode

        if (exitCode == 0):
            log('Executed command:\n{}\n'.format(command) +
                'Exit Code: {}\nOutput:\n{}'.format(exitCode, output))
        else:
            error_msg = 'Error running command:\n{}\n'.format(command)
            error_msg += 'Exit Code: {}\nOutput:\n{}'.format(exitCode, output)
            raise ValueError(error_msg)

    def genome_to_inputs(self, genome):
        """
        genome_to_inputs: convert genome annotation data (~json) to file inputs required by vConTACT
        :param genome:
        :return:
        """

        records = []
        gene2genome = OrderedDict()

        genome_data = genome['genomes'][0]

        for item in genome_data['data']['features']:
            if 'id' not in item:
                continue
                print('This feature does not have a valid id')
            elif 'dna_sequence' not in item or 'protein_translation' not in item:
                continue
                print('This feature {} does not have a valid DNA sequence.'.format(item['id']))
            else:
                # Create FASTA file
                if item['type'] == 'gene':
                    gene_record = SeqRecord(Seq(item['protein_translation'], IUPAC.protein), id=item['id'],
                                            description=item['function'])
                    records.append(gene_record)

                    # Build gene2genome
                    gene2genome.update({
                        item['id']: {
                            # 'contig_id': genome_data['data']['contig_ids'][0],
                            'contig_id': item['location'][0][0],
                            'protein_id': item['id'],
                            'keywords': item['function']
                        }
                    })

        return gene2genome, records

    def write_inputs(self, mapping, sequences):

        fasta_for_proteins_fp = os.path.join(self.scratch, 'vConTACT_proteins.fasta')
        with open(fasta_for_proteins_fp, 'w') as fasta_for_proteins_fh:
            SeqIO.write(sequences, fasta_for_proteins_fh, 'fasta')

        genes_to_genomes_mapping_fp = os.path.join(self.scratch, 'vConTACT_gene2genome.csv')
        with open(genes_to_genomes_mapping_fp, 'w') as genes_to_genomes_mapping_fh:
            fields = ['contig_id', 'protein_id', 'keywords']
            writer = csv.DictWriter(genes_to_genomes_mapping_fh, fieldnames=fields)
            writer.writeheader()

            for gene in mapping.keys():
                writer.writerow(mapping[gene])

        return genes_to_genomes_mapping_fp, fasta_for_proteins_fp

    def _generate_report(self, params):
        """
        _generate_report: generate summary report

        This will contain ALL the logic to generate the report, including areas that should/will be re-factored later

        """

        # Get filepath of summary file
        summary_fp = os.path.join(os.getcwd(), 'outdir', 'node_table_summary.csv')

        summary_df = pd.read_csv(summary_fp, header=0, index_col=0)
        html = summary_df.to_html(index=False, classes='my_class" id = "my_id')

        # Need to file write below
        direct_html = html_template.substitute(html_table=html)

        # html_dir = {
        #     'path': html_dir_path,
        #     'name': 'index.html',  # MUST match the filename of your main html page
        #     'description': 'My HTML report'
        # }

        report_params = {'message': 'Basic message to show in the report',
                         'report_object_name': 'import_matrix_from_excel_',
                         # Don't use until data objects that are created as result of running app
                         # 'objects_created': [{'ref': matrix_obj_ref,
                         #                      'description': 'Imported Matrix'}],
                         # Don't use until have files to attach to report
                         # 'file_links': [{}],
                         # Raw HTML
                         'direct_html': [direct_html],
                         'workspace_name': params['workspace_name'],
                         # 'html_links': [html_dir],
                         'direct_html_link_index': 0,
                         }

        kbase_report_client = KBaseReport(params['SDK_CALLBACK_URL'], token=params['KB_AUTH_TOKEN'])
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        # https://stackoverflow.com/a/600612/643675
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise
