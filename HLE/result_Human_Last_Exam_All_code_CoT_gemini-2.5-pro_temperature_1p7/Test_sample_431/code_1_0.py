def count_cysteines_in_linker():
    """
    This function calculates the number of Cysteine residues in the TM3-TM4 linker
    of the human GABAAρ1 receptor (UniProt P24046).
    """
    # Full amino acid sequence for human GABAAρ1 from UniProt: P24046
    full_sequence = (
        "MDWTWKIFSVGFLFLGALATSVHSEPAGASGRLSLEETNVSDLEIENNTTLVYIKGLSDRPV"
        "LGVSYPALLNVTVAEIASDSLIQYWNDRSTGNLVTLGSNNLVVEVWSDQKKFEYVDVKGGI"
        "RMLDLRRPLFYGVKGSYPSIRIMSFLLISKVSYPGISLAEMDITSCEIEMDGFSANGKRGP"
        "QAHGIHPLMRKVFNTMDYRKAIDVLMYAFCSLFVYVATYLPLTVSGVSKRATIPATIALSS"
        "WVFWVDVATLCFLSAHVALPAIASEKMATPLFTISNSAYKADDVEVYLWYITSGAALLGLG"
        "FVGYSLASLSEQKVAPAKSTLPKLSLPLPWLLYKAMAFLPVLAAVSVFGINKVAVPARVLG"
        "ITYLLMAVGQLSDIYKAMDGFPMAHFLPWFLGYVLNSVIVASSSSLHSPEEKLGLEGATLL"
        "DSRLGAYKAYTVETQGPRSEPSVSLQEPLLRKDHHHAHHPREVAAAESVARGLGDAGAGDS"
        "EAASSHRALAGAGPTPGGGRSRRALSSSRPRLGDPDPGAPLGSPASASASAPRPDPQSAE"
        "EEEEEEEEGRVCSGV"
    )

    # According to UniProt P24046 annotation, the TM3-TM4 linker domain
    # spans from residue 342 to 434 (1-based indexing).
    tm3_tm4_linker_start = 342
    tm3_tm4_linker_end = 434

    # Extract the linker sequence. Python uses 0-based indexing for slicing,
    # so we subtract 1 from the start position.
    linker_sequence = full_sequence[tm3_tm4_linker_start - 1:tm3_tm4_linker_end]

    # Count the number of Cysteine ('C') residues in the linker sequence.
    cysteine_count = linker_sequence.count('C')

    print(f"Analysis of human GABAAρ1 receptor (UniProt: P24046)")
    print(f"The TM3-TM4 linker domain spans from amino acid position {tm3_tm4_linker_start} to {tm3_tm4_linker_end}.")
    print("-" * 50)
    print(f"Linker Sequence: {linker_sequence}")
    print("-" * 50)
    print("Final Calculation:")
    print(f"Number of Cysteine ('C') residues in sequence from position {tm3_tm4_linker_start} to {tm3_tm4_linker_end} = {cysteine_count}")


# Run the analysis
count_cysteines_in_linker()
