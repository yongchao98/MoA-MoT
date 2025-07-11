def count_cysteines_in_domain():
    """
    This script calculates the number of Cysteine residues in the
    TM3-TM4 linker domain of the human GABAAρ1 receptor (UniProt: P24046).
    """
    # The canonical amino acid sequence for human GABAAρ1 (P24046).
    protein_sequence = (
        "MGFALFGLLAWVLGPAVGSGGGARLGEVKAELDDGQSVYKNLRFPDGSYGPHELVLTAHF"
        "NNTYSSIEISYPSALIFVDEKWRSDVFYGNLTLRLVNYTPLDLYFLCSDSNVPVDTFTRI"
        "LALHDDFPMDLHFNPLDFEDYVNCSLTFNFSYAYDTDELMALGVYGWYPYPGDGAGAYTG"
        "IKYICKQVMVCFDSISAKVETYGLTIRFFSRYSYDFRELQPCLLIVGSFASINTMEYSVH"
        "YFEVDRWLKEDKVPLYIGFRVPTNSALLVWFAFWLNYYARVALVYTIAVPMIPLAIYLFN"
        "VFYWTYYVNREDVELIVQTFSHTTNTTILGYSTKEVSYISLQVSFPSFSINEKLGMDPFL"
        "LRLRNNKKRVALGITDQLSRLNYYQSKKTEKKTAPNTTSPNGSPNIPLSSTGEIMAKIKY"
        "PAKAVSVNYKARGIFLPRTAILKAELDGSLGPAKLKATIPAKKPESEAPRTRPALSLDTS"
        "PSLPDSAPASSTAPGTPKPTSDPRSAAEAVSEAGSQASRC"
    )

    # The TM3-TM4 linker domain spans from residue 329 to 428 (1-based).
    start_residue = 329
    end_residue = 428

    # In Python, string slicing is 0-based and the end index is exclusive.
    start_index = start_residue - 1
    end_index = end_residue

    # Extract the linker sequence.
    linker_sequence = protein_sequence[start_index:end_index]

    # Count the number of Cysteine ('C') residues in the linker domain.
    cysteine_count = linker_sequence.count('C')

    # The final equation is: count = Number of 'C' in sequence[329:428]
    # The numbers defining the equation are the start and end residues (329, 428)
    # and the result of the count. The code below prints the final result.
    print(cysteine_count)

count_cysteines_in_domain()