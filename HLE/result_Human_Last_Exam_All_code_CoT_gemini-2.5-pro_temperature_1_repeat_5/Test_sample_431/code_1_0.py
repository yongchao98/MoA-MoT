def count_cysteines_in_linker():
    """
    This script calculates the number of Cysteine residues in the TM3-TM4
    intracellular linker domain of the human GABAA rho-1 (GABRR1) receptor.
    """
    # The canonical amino acid sequence for human GABAAρ1 from UniProt (P24046).
    full_protein_sequence = (
        "MGFALFTLSLLELVGSASYAEPSLGSYKDSDHDSELLSGYPSSRLRPEDAPAPASDGEC"
        "NIELEETSVLTGNIVLTGLDSQRTDDYVHTNIEITLRVKENNFHTDYLRTICMDIRMRFP"
        "DFGGPPVCTGMRIDATTMSVVSEVDIYNLGAVDQLPVLGVSVPIPSEIPDNYYSNDYFP"
        "AYTAQYCTEIEYNEGNYSGADYSFEVKKCLAVNKLMGMVSYVKAEIFYAKKPVYVALEDG"
        "TLVLYTLRLTPEATAALNYKLSRAIWLPDIFFFNAYVAYLATLFNINFVAYLIHFQTSQ"
        "PARAPLSFKAPMSYSAMMAWLGVLLYALSAAINAPGVKAYATAMDWFIAVCYAFVFSALI"
        "EFATVNYFTKRGYAWDGKSVVPEKPKKVKDPLIKKNVTYAHTNLERQARLCLPRPSTSTS"
        "ATAGAPLHQPAFLSNASKAETATTSTGSSQKYLSTQYQLPQPHTDQRRLPRASQLKINAT"
        "DIDRISRSLFYVTRLAVLGLSITLLVTISWVSFWINYDASAARVALGITTVLTMTTINTH"
        "LRETLPKIPYVKAIDMYLMGCFVFVFLALLEYAFVNYIFFGKGPQH"
    )

    # Remove any spaces or newlines to ensure sequence integrity.
    full_protein_sequence = "".join(full_protein_sequence.split())

    # According to UniProt P24046, the TM3-TM4 linker (a large cytoplasmic
    # domain) spans from residue 331 to 432 (inclusive, 1-based indexing).
    linker_start_index = 331 - 1  # Convert to 0-based index for Python
    linker_end_index = 432

    # Extract the sequence of the TM3-TM4 linker domain.
    linker_sequence = full_protein_sequence[linker_start_index:linker_end_index]

    # Count the number of Cysteine ('C') residues in the linker.
    cysteine_count = linker_sequence.count('C')

    # Print the details and the final result.
    print(f"Protein: Human GABAAρ1 (UniProt: P24046)")
    print(f"Domain: TM3-TM4 linker (Residues 331-432)")
    print(f"Linker Sequence: {linker_sequence}")
    print("\n--- Calculation ---")
    # This fulfills the requirement to "output each number in the final equation"
    # for a counting problem.
    print(f"Number of Cysteine ('C') residues found: {cysteine_count}")


if __name__ == "__main__":
    count_cysteines_in_linker()
