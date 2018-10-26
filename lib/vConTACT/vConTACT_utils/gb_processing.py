from Bio import SeqIO
import pandas as pd
from io import StringIO

gbk_txt = """
##sequence-region NC_003287.2
NC_003287.2	KBase	gene	6006	6407	.	+	0	ID=M13p01; db_xref=GeneID:927328; locus_tag=M13p01; protein_id=NP_510885.1; gene=II; product=replication initiation protein
NC_003287.2	KBase	gene	1	831	.	+	0	ID=M13p01; db_xref=GeneID:927328; locus_tag=M13p01; protein_id=NP_510885.1; gene=II; product=replication initiation protein
NC_003287.2	KBase	CDS	6006	6407	.	+	0	ID=M13p01_CDS_1; Parent=M13p01; note=The pII protein binds to the intergenic region in double-stranded%2C supercoiled DNA %28replicative form%29 and nicks the positive strand. The resulting free 3%27 OH terminus serves as a primer for new strand synthesis using host factors %28Rep helicase%2C DNA polymerase III%29. After a full round of replication%2C pII circularizes the displaced DNA strand. Newly synthesized %28%2B%29 single stranded DNA is either converted to RF DNA or packaged as genomes%2C depending on the concentrations of ssDNA binding protein pV and a separately expressed fragment of pII termed pX. Functions: DNA replication initiation; db_xref=GOA:P03659; db_xref=GeneID:927328; db_xref=UniProtKB/Swiss-Prot:P03659; gene=II; locus_tag=M13p01; protein_id=NP_510885.1; product=replication initiation protein
NC_003287.2	KBase	CDS	1	831	.	+	0	ID=M13p01_CDS_1; Parent=M13p01; note=The pII protein binds to the intergenic region in double-stranded%2C supercoiled DNA %28replicative form%29 and nicks the positive strand. The resulting free 3%27 OH terminus serves as a primer for new strand synthesis using host factors %28Rep helicase%2C DNA polymerase III%29. After a full round of replication%2C pII circularizes the displaced DNA strand. Newly synthesized %28%2B%29 single stranded DNA is either converted to RF DNA or packaged as genomes%2C depending on the concentrations of ssDNA binding protein pV and a separately expressed fragment of pII termed pX. Functions: DNA replication initiation; db_xref=GOA:P03659; db_xref=GeneID:927328; db_xref=UniProtKB/Swiss-Prot:P03659; gene=II; locus_tag=M13p01; protein_id=NP_510885.1; product=replication initiation protein
NC_003287.2	KBase	gene	496	831	.	+	0	ID=M13p02; db_xref=GeneID:927329; locus_tag=M13p02; gene=X; protein_id=NP_510886.1; product=gene X product
NC_003287.2	KBase	CDS	496	831	.	+	0	ID=M13p02_CDS_1; Parent=M13p02; note=Protein pX is synthesized from an in-frame internal initiation codon in gene II and is identical to the C-terminal third of pII. Required for the accumulation of %28%2B%29 single stranded genomic DNA.; db_xref=GOA:P03659; db_xref=GeneID:927329; db_xref=UniProtKB/Swiss-Prot:P03659; gene=X; locus_tag=M13p02; protein_id=NP_510886.1; product=gene X product
NC_003287.2	KBase	gene	843	1106	.	+	0	ID=M13p03; db_xref=GeneID:927330; protein_id=NP_510887.1; locus_tag=M13p03; gene=V; product=helix destabilising protein
NC_003287.2	KBase	CDS	843	1106	.	+	0	ID=M13p03_CDS_1; Parent=M13p03; db_xref=GOA:P03669; db_xref=GeneID:927330; db_xref=UniProtKB/Swiss-Prot:P03669; gene=V; locus_tag=M13p03; protein_id=NP_510887.1; product=helix destabilising protein
NC_003287.2	KBase	old_sequence	898	898	.	+	0	ID=M13p03_old_sequence_1; Parent=M13p03; gene=V; locus_tag=M13p03"""


def gbk_to_inputs(gbk_fp):
    """

    :param gbk_fp:
    :return:
    """

    with open(gbk_fp, 'rU') as gbk_fh:

        for record in SeqIO.parse(gbk_fh, 'genbank'):
            print(dir(record))


# gbk_to_inputs('/Users/bolduc.10/EMBOSS-6.6.0/test/data/seqexample.gff3')
gbk_to_inputs('/Users/bolduc.10/EMBOSS-6.6.0/test/data/featexample2.gff3')

"""
or record in SeqIO.parse(genbank_fh, 'genbank'):
        for f in record.features:
            if f.type == 'CDS':
                newRecord = SeqRecord(Seq(f.qualifiers['translation'][0], IUPAC.protein),
                                      id=f.qualifiers['locus_tag'][0],
                                      description=f.qualifiers['product'][0])"""