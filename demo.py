import streamlit as st
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio import pairwise2
from Bio.Seq import Seq
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt
import os

# Title of the Streamlit App
st.title("🔬 Phylogenetic Tree & Sequence Analysis")

# Input box for FASTA sequences
st.subheader("Enter FASTA Sequences")
fasta_input = st.text_area(
    "Paste sequences in FASTA format:",
    """>Seq1
ATGCGTACGTTAG
>Seq2
ATGCGTACGTGAG
>Seq3
ATGCGTTCGTTAG
>Seq4
ATGCGGACGTTAG""",
    height=200,
)


# Function to save input sequences to a FASTA file
def save_fasta(input_text, file_name="sequences.fasta"):
    with open(file_name, "w") as f:
        f.write(input_text.strip())


# Function to analyze sequences (alignment, SNPs, identity)
def analyze_sequences(fasta_file, reference_seq):
    alignment_results = []
    alignment = AlignIO.read(fasta_file, "fasta")

    for record in alignment:
        seq = str(record.seq)
        alignments = pairwise2.align.globalxx(reference_seq, seq)
        best_alignment = alignments[0]
        similarity = sum(a == b for a, b in zip(best_alignment.seqA, best_alignment.seqB))
        percent_identity = (similarity / len(best_alignment.seqA)) * 100

        # Convert SNPs list of tuples to a string for better handling in DataFrame
        snps = [(i + 1, ref, seq[i]) for i, ref in enumerate(reference_seq) if seq[i] != ref]
        snps_str = "; ".join([f"Pos {pos}: {ref} → {alt}" for pos, ref, alt in snps])

        alignment_results.append({
            "Sequence": record.id,
            "Percent Identity": percent_identity,
            "SNPs": snps_str,  # Store SNPs as a string instead of a list
            "Best Alignment": f"{best_alignment.seqA}\n{best_alignment.seqB}"
        })

    return pd.DataFrame(alignment_results)


# Function to build phylogenetic tree
def build_phylogenetic_tree(fasta_file):
    alignment = AlignIO.read(fasta_file, "fasta")
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, "upgma")
    tree = constructor.build_tree(alignment)

    # Save tree as an image
    fig = plt.figure(figsize=(6, 4), dpi=100)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax)
    plt.savefig("tree.png")
    return "tree.png"


# Run analysis when the button is clicked
if st.button("Generate Phylogenetic Tree & Analyze Sequences"):
    if not fasta_input.strip():
        st.error("Please enter sequences in FASTA format!")
    else:
        # Save user input to FASTA file
        fasta_file = "sequences.fasta"
        save_fasta(fasta_input, fasta_file)
        st.success("FASTA file saved successfully!")

        # Reference sequence (first sequence in input)
        first_seq = fasta_input.split("\n")[1].strip()

        # Generate phylogenetic tree
        st.subheader("📌 Phylogenetic Tree")
        try:
            tree_image = build_phylogenetic_tree(fasta_file)
            st.image(tree_image, caption="Phylogenetic Tree")
        except Exception as e:
            st.error(f"Error generating tree: {e}")

        # Analyze sequences
        st.subheader("📊 Sequence Analysis")
        try:
            sequence_results = analyze_sequences(fasta_file, first_seq)
            st.dataframe(sequence_results)
        except Exception as e:
            st.error(f"Error analyzing sequences: {e}")
