import time
import os
import subprocess
import csv
import uuid
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from string import Template
from pyparsing import Literal, SkipTo
import pandas as pd

from KBaseReport.KBaseReportClient import KBaseReport
from Workspace.WorkspaceClient import Workspace
from GenomeAnnotationAPI.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from DataFileUtil.DataFileUtilClient import DataFileUtil as dfu
from installed_clients.AssemblyUtilClient import AssemblyUtil
# from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil as gfu
from KBaseDataObjectToFileUtils.KBaseDataObjectToFileUtilsClient import KBaseDataObjectToFileUtils as ofu


def log(message, prefix_newline=False):
    """
    Logging function, provides a hook to suppress or redirect log messages.
    """
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


html_template = Template("""<!DOCTYPE html>
<html lang="en">
  <head>
  
    <link href="https://netdna.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdn.datatables.net/1.10.22/css/jquery.dataTables.min.css" rel="stylesheet">
    <link href="https://cdn.datatables.net/buttons/1.5.2/css/buttons.dataTables.min.css" rel="stylesheet">

    <link href="https://cdn.datatables.net/searchpanes/1.2.0/css/searchPanes.dataTables.min.css" rel="stylesheet">
    <link href="https://cdn.datatables.net/select/1.3.1/css/select.dataTables.min.css" rel="stylesheet">

    <script src="https://code.jquery.com/jquery-3.5.1.js" type="text/javascript"></script>
    <script src="https://cdn.datatables.net/1.10.22/js/jquery.dataTables.min.js" type="text/javascript"></script>
    <script src="https://cdn.datatables.net/buttons/1.6.4/js/dataTables.buttons.min.js" type="text/javascript"></script>
    <script src="https://cdn.datatables.net/buttons/1.6.4/js/buttons.flash.min.js" type="text/javascript"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js" type="text/javascript"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/pdfmake.min.js" type="text/javascript"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/vfs_fonts.js" type="text/javascript"></script>
    <script src="https://cdn.datatables.net/buttons/1.6.4/js/buttons.html5.min.js" type="text/javascript"></script>
    <script src="https://cdn.datatables.net/buttons/1.6.4/js/buttons.print.min.js" type="text/javascript"></script>

    <script src="https://cdn.datatables.net/searchpanes/1.2.0/js/dataTables.searchPanes.min.js" type="text/javascript"></script>
    <script src="https://cdn.datatables.net/select/1.3.1/js/dataTables.select.min.js" type="text/javascript"></script>
  
    <style>
    tfoot input {
        width: 100%;
        padding: 3px;
        box-sizing: border-box;
    }
    </style>
  
  </head>
  
  <body>
  
    <div class="container">
      <div>
        ${html_table}
      </div>
    </div>

    <script type="text/javascript">
      $$(document).ready(function() {
        $$('#my_id tfoot th').each( function () {
          var title = $$(this).text();
          $$(this).html( '<input type="text" placeholder="Search '+title+'" />' );
        });

        var table = $$('#my_id').DataTable({
          buttons: [
            'copy', 'csv', 'excel', 'pdf', 'print'],
          scrollX: true,
          dom: 'lPfrtip'  //Necessary for buttons to work
        });

        table.columns().every( function () {
          var that = this;

          $$( 'input', this.footer() ).on( 'keyup change', function () {
            if ( that.search() !== this.value ) {
              that
              .search( this.value )
              .draw();
            }
          });
        } );
      } );
    </script>
    
  </body>
</html>""")


class vConTACTUtils:

    def __init__(self, config):
        self.scratch = os.path.abspath(config['scratch'])
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.token = os.environ['KB_AUTH_TOKEN']
        self.scratch = os.path.abspath(config['scratch'])
        self.ws = Workspace(config['workspace-url'], token=self.token)
        self.genome_api = GenomeAnnotationAPI(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)

    def vcontact_help(self):
        command = "vcontact --help"
        self._run_command(command)

    def execute(self, command: list):
        """
        :param command: Command suitable for running in subprocess, must use a ['ls', '-l'] format
        :return: Response from command
        """
        # logger.info('Running command: {}'.format(command))
        print('Running command: {}'.format(' '.join(command)))
        res = subprocess.run(command, shell=False, encoding='utf-8', check=True)

        return res

    def run_vcontact(self, params):

        # Determine KBase "inputs" for vConTACT2
        genome = params['genome']

        obj_type = self.ws.get_object_info3({'objects': [{'ref': genome}]})['infos'][0][2]

        if 'assembly' in obj_type.lower():  # If KBaseGenomeAnnotations.Assembly

            # Assembly requires annotation
            genome_fp = self.au.get_assembly_as_fasta({'ref': genome})['path']
            proteins_fp = os.path.join(self.scratch, 'proteins.faa')
            proteins_gbk = os.path.join(self.scratch, 'proteins.gbk')
            gene2genome_fp = os.path.join(self.scratch, 'gene2genome.csv')

            prodigal_cmd = ['prodigal', '-a', proteins_fp, '-o', proteins_gbk, '-f', 'gbk',
                            '-i', genome_fp, '-p', 'meta']
            res = self.execute(prodigal_cmd)

            records = {}
            with open(proteins_fp, 'r') as proteins_fh:
                for record in SeqIO.parse(proteins_fh, 'fasta'):

                    records[len(records)] = {
                        'protein_id': record.id,
                        'contig_id': record.id.rsplit('_', 1)[0],
                        'keywords': 'None'
                    }

            g2g_df = pd.DataFrame.from_dict(records, orient='index')
            g2g_df.to_csv(gene2genome_fp, index=False)

            # Pass filepaths to the app and run
            params['gene2genome'] = gene2genome_fp
            params['sequences'] = proteins_fp

        elif 'kbasegenomes' in obj_type.lower(): # If KBaseGenomes.Genome
            genome_data = self.genome_api.get_genome_v1({"genomes": [{"ref": genome}]})

            # Convert genome data into "reasonable" parse form and write to scratch filesystem
            gene2genome, sequences = self.genome_to_inputs(genome_data)
            gene2genome_fp, sequences_fp = self.write_inputs(gene2genome, sequences)

            # Pass filepaths to the app and run
            params['gene2genome'] = gene2genome_fp
            params['sequences'] = sequences_fp

        elif 'binnedcontigs' in obj_type.lower():  # If KBaseMetagenomes.BinnedContigs
            print('KBaseMetagenomes.BinnedContigs hasnt been enabled. Check back later.')
            exit(1)
        else:
            print('Unknown error in identifying object types')

        print('Available database files')
        print(os.listdir('/miniconda/lib/python3.7/site-packages/vcontact2/data/'))

        # Just iterate through all parameters
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
        command = 'vcontact2 --output-dir outdir'
        # Binaries
        command += ' --diamond-bin /usr/local/bin/diamond --c1-bin /usr/local/bin/cluster_one-1.0.jar'

        for param, cmd in mappings.items():
            command += ' {} {}'.format(cmd, params[param])

        self._run_command(command)

        report = self._generate_report(params)

        return report

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
                    desc = (item['functions'] if item.get('functions', None)
                                              else item.get('function', ''))
                    gene_record = SeqRecord(Seq(item['protein_translation']), id=item['id'],
                                            description=desc)
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

        # Get
        self.dfu = dfu(self.callback_url)

        # Get filepath of summary file
        summary_fp = os.path.join(os.getcwd(), 'outdir', 'genome_by_genome_overview.csv')

        summary_df = pd.read_csv(summary_fp, header=0, index_col=0)
        html = summary_df.to_html(index=False, classes='my_class table-striped" id = "my_id')

        # Need to file write below
        direct_html = html_template.substitute(html_table=html)

        # Find header so it can be copied to footer, as dataframe.to_html doesn't include footer
        start_header = Literal("<thead>")
        end_header = Literal("</thead>")

        text = start_header + SkipTo(end_header)

        new_text = ''
        for data, start_pos, end_pos in text.scanString(direct_html):
            new_text = ''.join(data).replace(' style="text-align: right;"', '').replace('thead>',
                                                                                        'tfoot>\n  ') + '\n</tfoot>'

        # Get start and end positions to insert new text
        end_tbody = Literal("</tbody>")
        end_table = Literal("</table>")

        insertion_pos = end_tbody + SkipTo(end_table)

        final_html = ''
        for data, start_pos, end_pos in insertion_pos.scanString(direct_html):
            final_html = direct_html[:start_pos + 8] + '\n' + new_text + direct_html[start_pos + 8:]

        output_dir = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_dir)
        result_fp = os.path.join(output_dir, 'index.html')

        with open(result_fp, 'w') as result_fh:
            result_fh.write(final_html)

        report_shock_id = self.dfu.file_to_shock({
            'file_path': output_dir,
            'pack': 'zip'
        })['shock_id']

        html_report = [{
            'shock_id': report_shock_id,
            'name': os.path.basename(result_fp),
            'label': os.path.basename(result_fp),
            'description': 'HTML summary report for vConTACT2'
        }]

        report_params = {'message': 'Basic message to show in the report',
                         'workspace_name': params['workspace_name'],
                         'html_links': html_report,
                         'direct_html_link_index': 0,
                         'report_object_name': 'vConTACT_report_{}'.format(str(uuid.uuid4())),
                         # Don't use until have files to attach to report
                         # 'file_links': [{}],
                         # Don't use until data objects that are created as result of running app
                         # 'objects_created': [{'ref': matrix_obj_ref,
                         #                      'description': 'Imported Matrix'}],
                         }

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
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
