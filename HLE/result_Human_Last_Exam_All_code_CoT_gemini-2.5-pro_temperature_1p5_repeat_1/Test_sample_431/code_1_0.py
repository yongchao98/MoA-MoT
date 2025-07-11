def count_cysteines_in_linker():
    """
    This script calculates the number of Cysteine residues in the TM3-TM4 linker
    of the human GABAA receptor subunit rho-1 (GABRR1).
    """

    # Full amino acid sequence for human GABRR1 from UniProt (P24046)
    full_sequence = (
        "MHFRLWGAWALALLSVPCLALAAPSSSAKIDRLYVGYDGKLRPDIGVKVMINNTLVMTTS"
        "IHARNSSLPDVSLYYNDMGFYAWAGGWDPEDDATNVTRSPLKRIAYTYPNTLVVSCVSPL"
        "HLEDFPMDAHACPLKFGSYAYPKDIWMWYVDASAIDLKKFVSYAYTTNCMFIVSSILYFA"
        "YAFNFNEAPVAKVTAPIGTYATMDWFFANGKKSVEHDITVFENKLTLDHKYAYTEIILHY"
        "WCCCYIVMSAIWMAFWINSRVPARVSLGITTVLTMTTLSTIARKSLPKVSYVKAIDIWMAV"
        "CLLAFVFSALIEFATVNYFTKRGYAWDGKSVVPEKPKKVKDPLIKKNMTILRAFPLLFGA"
        "LVLNSYYSRAYADQLKVSEVPTVALGVKTVLMRTSFPRLGVELGCTDIHICSYMGYFTIN"
        "TYMPCTLIVYNALVTGVSPRLSIHFRLRRNKPALLRISAPAKAKAEVPRSRSASLRERRL"
        "HSRRELKEYARAAALKSRAPTPHQAAALQVSWIGYVMLIAMVSSISAWLFWYEASAYARV"
        "ALGVYSAAIRMMLTLL"
    )

    # According to UniProt P24046, the TM3-TM4 linker (a large intracellular loop)
    # spans from residue 341 to 447.
    start_pos = 341
    end_pos = 447

    # Python slicing is 0-indexed, so we subtract 1 from the start position.
    # The end position in the slice is exclusive, so `end_pos` correctly extracts up to residue 447.
    linker_sequence = full_sequence[start_pos - 1:end_pos]

    # Count the number of Cysteine ('C') residues in the extracted linker sequence.
    cysteine_count = linker_sequence.count('C')

    print(f"Protein: Human GABAAρ1 receptor (UniProt P24046)")
    print(f"Domain of interest: TM3-TM4 linker")
    print(f"Residue range for domain: {start_pos}-{end_pos}")
    print(f"\nLinker Sequence ({len(linker_sequence)} residues):")
    print(linker_sequence)
    
    # "Output each number in the final equation" - showing the components of the count
    print(f"\nCounting 'C' in the sequence above...")
    print(f"Total number of Cysteine residues = {cysteine_count}")

# Execute the function
count_cysteines_in_linker()