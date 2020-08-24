import matplotlib.pyplot as plt
import biotite
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.io.fasta as fasta
import biotite.database.entrez as entrez
import biotite.sequence.graphics as graphics

# Download and parse protein sequences of Covid and Mers
covid_file_path = entrez.fetch("NC_045512","myresult_dir",suffix="fa",db_name="nuccore",ret_type="fasta")
mers_file_path = entrez.fetch("NC_019843.3","myresult_dir",suffix="fa",db_name="nuccore",ret_type="fasta")
# Read the file
c_file = fasta.FastaFile()
c_file.read(covid_file_path)
m_file = fasta.FastaFile()
m_file.read(mers_file_path)
# Display
for h,s in c_file.items():
    print(h)
    print(s)
    covid_seq = seq.NucleotideSequence(s)
for h,s in m_file.items():
    print(h)
    print(s)
    mers_seq = seq.NucleotideSequence(s)
mini_covid_seq = covid_seq[0:100]
mini_mers_seq = mers_seq[0:100]
matrix = align.SubstitutionMatrix.std_nucleotide_matrix()
# Perform pairwise sequence alignment with affine gap penalty
# Terminal gaps are not penalized
alignments = align.align_optimal(mini_covid_seq, mini_mers_seq, matrix,gap_penalty=(-10, -1), terminal_penalty=False)
                                
# Draw first and only alignment
# The color intensity indicates the similiarity
fig = plt.figure(figsize=(8.0, 2.5))
ax = fig.add_subplot(111)
graphics.plot_alignment_similarity_based(
    ax, alignments[0], matrix=matrix, labels=["SARS_Covid", "MERS"],
    show_numbers=True, show_line_position=True
)
fig.tight_layout()

plt.show()

# Draw first and only alignment
# The color intensity indicates the similiarity
fig = plt.figure(figsize=(8.0, 2.5))
ax = fig.add_subplot(111)
graphics.plot_alignment_similarity_based(
    ax, alignments[0], matrix=matrix, labels=["SARS_Covid", "MERS"],
    show_numbers=True, show_line_position=True
)
fig.tight_layout()

plt.show()
