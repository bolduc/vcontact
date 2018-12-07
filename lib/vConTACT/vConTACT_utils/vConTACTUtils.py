import time
import os
import subprocess
import csv
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def log(message, prefix_newline=False):
    """
    Logging function, provides a hook to suppress or redirect log messages.
    """
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


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
                print('This feature does not have a valid id')
            elif 'dna_sequence' not in item or 'protein_translation' not in item:
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
                            'contig_id': genome_data['data']['location'][0][0],
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
