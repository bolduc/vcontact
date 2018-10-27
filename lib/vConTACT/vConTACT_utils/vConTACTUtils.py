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

        print(type(genome))

        print(genome.keys())

        for item in genome['data']['features']:
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
                            'contig_id': genome['data']['contig_ids'][0],
                            'protein_id': item['id'],
                            'keywords': item['function']
                        }
                    })

        return gene2genome, records


    def write_inputs(self, mapping, sequences):

        fasta_for_proteins_fp = os.path.join(self.scratch, 'vConTACT_proteins.fasta')
        with open(fasta_for_proteins_fp, 'w') as fasta_for_proteins_fh:
            SeqIO.write(sequences, fasta_for_proteins_fh, 'fasta')

        genes_to_genomes_mapping_fp = os.path.join(self.scratch, 'vConTACT_gene2genome.fasta')
        with open(genes_to_genomes_mapping_fp, 'w') as genes_to_genomes_mapping_fh:
            fields = ['contig_id', 'protein_id', 'keywords']
            writer = csv.DictWriter(genes_to_genomes_mapping_fh, fieldnames=fields)
            writer.writeheader()

            for gene in mapping.keys():
                writer.writerow(mapping[gene])

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
