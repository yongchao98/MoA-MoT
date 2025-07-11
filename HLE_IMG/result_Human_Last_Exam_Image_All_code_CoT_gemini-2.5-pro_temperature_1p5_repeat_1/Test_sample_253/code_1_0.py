def count_amino_acids():
    """
    This function calculates the number of amino acids in the XPH1 protein
    from the lichen Xanthoria parietina.
    """
    # The amino acid sequence for Xanthoria parietina hydrophobin 1 (XPH1)
    # was retrieved from the NCBI database (Accession: CAC84045.1).
    protein_sequence = "MALLTSLLAALAVAPVLAATSPTTCTDGFNGCKPFCCAADKARGKLLACCSTSLPCTNPGPYHTHCGTSGCIPHACCPTIGSMLSCCPVPAGL"

    # Each letter in the sequence represents one amino acid.
    # We calculate the total number by finding the length of the string.
    number_of_amino_acids = len(protein_sequence)

    print(f"The XPH1 protein of Xanthoria parietina has the following amino acid sequence:\n{protein_sequence}\n")
    print(f"To find the total number of amino acids, we calculate the length of this sequence.")
    print(f"The length of the XPH1 protein is {number_of_amino_acids} amino acids.")

# Run the function to print the result
count_amino_acids()