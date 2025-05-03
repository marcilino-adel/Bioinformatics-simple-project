from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np


class ORFFinder:
    """Class to detect Open Reading Frames (ORFs) in a DNA sequence."""

    @staticmethod
    def find_orfs(sequence: Seq) -> list:
        """
        Find all ORFs in a DNA sequence.

        Args:
            sequence: DNA sequence (Bio.Seq.Seq).

        Returns:
            List of tuples (orf, start, end) for each ORF, where:
                - orf: ORF sequence (str).
                - start: 0-based start position (int).
                - end: 0-based end position (int, exclusive).
        """
        # Initialize list to store ORFs
        orfs = []
        # Convert sequence to string for processing
        seq_str = str(sequence)
        # Iterate over three reading frames
        for frame in range(3):
            # Check each codon position in the frame
            for i in range(frame, len(seq_str) - 2, 3):
                # Look for start codon (ATG)
                if seq_str[i:i + 3] == "ATG":
                    # Search for stop codon
                    for j in range(i, len(seq_str) - 2, 3):
                        codon = seq_str[j:j + 3]
                        # If stop codon found, extract ORF
                        if codon in ["TAA", "TAG", "TGA"]:
                            orf = seq_str[i:j + 3]
                            orfs.append((orf, i, j + 3))
                            break
        return orfs

    @staticmethod
    def find_orfs_modified(sequence: Seq, min_length: int) -> list:
        """
        Find ORFs in a DNA sequence, filtering out those shorter than min_length.

        Args:
            sequence: DNA sequence (Bio.Seq.Seq).
            min_length: Minimum ORF length in nucleotides (int).

        Returns:
            List of tuples (orf, start, end) for ORFs meeting length requirement.
        """
        # Initialize list to store ORFs
        orfs = []
        # Convert sequence to string for processing
        seq_str = str(sequence)
        # Iterate over three reading frames
        for frame in range(3):
            # Check each codon position in the frame
            for i in range(frame, len(seq_str) - 2, 3):
                # Look for start codon (ATG)
                if seq_str[i:i + 3] == "ATG":
                    # Search for stop codon
                    for j in range(i, len(seq_str) - 2, 3):
                        codon = seq_str[j:j + 3]
                        # If stop codon found, check ORF length
                        if codon in ["TAA", "TAG", "TGA"]:
                            orf = seq_str[i:j + 3]
                            # Only include ORFs meeting minimum length
                            if len(orf) >= min_length:
                                orfs.append((orf, i, j + 3))
                            break
        return orfs


class DotPlotGenerator:
    """Class to generate dot plots for comparing two DNA sequences."""

    @staticmethod
    def dot_plot(seq1: Seq, seq2: Seq, output_file: str = "dot_plot.png") -> None:
        """
        Generate a dot plot for two DNA sequences, placing a dot where nucleotides match.
        Saves the plot to a file and displays it.

        Args:
            seq1: First DNA sequence (Bio.Seq.Seq).
            seq2: Second DNA sequence (Bio.Seq.Seq).
            output_file: File to save the plot (default: 'dot_plot.png').
        """
        # Initialize lists to store matching positions
        x, y = [], []
        # Compare each nucleotide pair
        for i, base1 in enumerate(seq1):
            for j, base2 in enumerate(seq2):
                # If nucleotides match, record position
                if base1 == base2:
                    x.append(i + 0.5)  # Offset for better visualization
                    y.append(j + 0.5)

        # Create scatter plot
        plt.figure(figsize=(8, 8))
        plt.scatter(x, y, s=300, c='black', marker='*')
        plt.title("Dot Plot of Two DNA Sequences")
        plt.xlabel("Sequence 1")
        plt.ylabel("Sequence 2")
        plt.grid(True)
        plt.savefig(output_file)
        plt.show()
        plt.close()


class CodonUsageAnalyzer:
    """Class to analyze codon usage frequencies in DNA sequences."""

    @staticmethod
    def codon_usage_table(sequence: Seq) -> dict:
        """
        Generate a codon usage frequency table for a single DNA sequence.
        Returns a dictionary of codon counts and saves a bar plot.

        Args:
            sequence: DNA sequence (Bio.Seq.Seq).

        Returns:
            Dictionary with codons as keys and their counts as values.
        """
        # Initialize codon counter
        codons = Counter()

        # Convert sequence to string and trim to be divisible by 3
        seq_str = str(sequence)
        if len(seq_str) % 3 != 0:
            seq_str = seq_str[:-(len(seq_str) % 3)]

        # Count codons in-frame
        for i in range(0, len(seq_str) - 2, 3):
            codon = seq_str[i:i + 3]
            if len(codon) == 3:
                codons[codon] += 1

        # Visualize as bar plot
        plt.figure(figsize=(12, 6))
        plt.bar(codons.keys(), codons.values())
        plt.title("Codon Usage Frequency")
        plt.xlabel("Codon")
        plt.ylabel("Count")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig("codon_usage.png")
        plt.show()
        plt.close()

        return dict(codons)

    @staticmethod
    def codon_usage_table_modified(sequences: list, output_file: str = "codon_usage_combined.png") -> dict:
        """
        Generate a combined codon usage frequency table for multiple DNA sequences.
        Returns a dictionary of average codon frequencies and saves a bar plot.

        Args:
            sequences: List of DNA sequences (Bio.Seq.Seq).
            output_file: File to save the plot (default: 'codon_usage_combined.png').

        Returns:
            Dictionary with codons as keys and their average frequencies as values.
        """
        # Initialize combined codon counter
        all_codons = Counter()
        total_codons = 0

        # Process each sequence
        for seq in sequences:
            seq_str = str(seq)
            if len(seq_str) % 3 != 0:
                seq_str = seq_str[:-(len(seq_str) % 3)]
            for i in range(0, len(seq_str) - 2, 3):
                codon = seq_str[i:i + 3]
                if len(codon) == 3:
                    all_codons[codon] += 1
                    total_codons += 1

        # Calculate average frequencies (normalize by total codons)
        codon_freq = {codon: count / total_codons for codon, count in all_codons.items()}

        # Visualize as bar plot
        plt.figure(figsize=(12, 6))
        plt.bar(codon_freq.keys(), codon_freq.values())
        plt.title("Combined Codon Usage Frequency (Normalized)")
        plt.xlabel("Codon")
        plt.ylabel("Frequency")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_file)
        plt.show()
        plt.close()

        return codon_freq


def main():
    """
    Main function to demonstrate usage of ORF detection, dot plot, and codon usage analysis.
    Includes example usage with string sequences and FASTA file input.
    """
    # Point 1: ORF Detection
    try:
        # Read the first sequence from FASTA file
        records = SeqIO.parse("exFasta.fasta", "fasta")
        dna = next(records).seq
        # Find ORFs with minimum length
        orfs = ORFFinder.find_orfs_modified(dna, min_length=30)
        print("ORFs found:")
        for orf, start, end in orfs:
            print(f"ORF: {orf}, Start: {start}, End: {end}")
    except StopIteration:
        print("Error: The FASTA file is empty (no sequences found).")
    except FileNotFoundError:
        print("Error: The file 'exFasta.fasta' was not found.")

    # Point 4: Dot Plot
    seq1 = Seq("ATCGATCG")
    seq2 = Seq("ATAGCTCG")
    DotPlotGenerator.dot_plot(seq1, seq2, "dot_plot_string.png")
    print("Dot plot saved as: dot_plot_string.png")

    # Point 5: Codon Usage (Single Sequence)
    seq = Seq("ATGCGTAAATAGATGCGT")
    codon_freq = CodonUsageAnalyzer.codon_usage_table(seq)
    print("\nCodon Usage Table (Single Sequence):")
    for codon, count in codon_freq.items():
        print(f"{codon}: {count}")
    print("Codon usage plot saved as: codon_usage.png")

    # Point 5: Codon Usage (Multiple Sequences from FASTA)
    try:
        seqs = [record.seq for record in SeqIO.parse("exFasta.fasta", "fasta")]
        codon_freq = CodonUsageAnalyzer.codon_usage_table_modified(seqs, "codon_usage_combined_fasta.png")
        print("\nCombined Codon Usage Table (FASTA File):")
        for codon, freq in codon_freq.items():
            print(f"{codon}: {freq:.4f}")
        print("Combined codon usage plot saved as: codon_usage_combined_fasta.png")
    except FileNotFoundError:
        print("FASTA file not found. Please provide exFasta.fasta.")


if __name__ == "__main__":
    main()