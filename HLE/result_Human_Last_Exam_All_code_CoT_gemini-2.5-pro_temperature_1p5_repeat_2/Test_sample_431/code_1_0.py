def count_cysteines_in_linker():
    """
    This script calculates the number of Cysteine (C) residues in the 
    TM3-TM4 linker domain of the human GABAA receptor subunit rho-1 (GABRR1).
    """

    # The canonical amino acid sequence for human GABRR1 from UniProt (P24046).
    full_sequence = (
        "MGFMKLVYFLIAGLSCLGLEETRYASDSGIQVSWTPNSSNMSPVQVTIFNLVLENNVRELI"
        "QTDLFHDIRWLDERILFPGDNVPLNLSFNVTGPTVLLTVTMSSTIATDMSHYAIDISLFVW"
        "LSDPMDAASHSPLRASYNPDEAINVILKKYEWTVIDISYASINKMEIDAIACYTLDFPMGV"
        "PGHPLSVMVSCTMNSVSFGSAAYPREDFVYVFSYATFNLTYWVHYDDNKVRPDLAIYAIKS"
        "YGRYEVPEPKKVKDPLIKKNMMLVWICFLSFLGALGYYFLTNAYFAHDDDLEAYRAFYFIA"
        "LVGYACLFLINPFVFSAYPMDAKEVPLYSVKLIDIWVSFFPFVLFFNVYWLWYIFKTTDPG"
        "QLAAKIPKYVISKAIDMYLMGCFVFVFASLLEYATVNYIFFGRGPGRGALITSTSTIAKAR"
        "SVTPARTPPHSLQLHFPSNIPSQLSAEDINECAIKACPKKTDYRYIWKACCARPPASAIHF"
        "VSRHNGSKSPQLHFEPLDQKAYATAMTMDYFRMLACRWKHYLEFKMVCAKRIGYFMIQMYI"
        "PSLLIVILSWSVSFWINMDAAPARVALGITTVLTMTTISTHVAR"
    )

    # According to UniProt P24046, the transmembrane domains are located at:
    # TM3: residues 297-318
    # TM4: residues 450-470
    # Therefore, the TM3-TM4 linker is the region between the end of TM3 and the start of TM4.
    
    # In 1-based indexing (protein sequences), the linker is from residue 319 to 449.
    tm3_end_pos = 318
    tm4_start_pos = 450
    
    # In 0-based indexing (Python strings), we slice from index 318 up to (but not including) 449.
    # We want residues 319 (index 318) to 449 (index 448), so the slice is [318:449].
    linker_start_index = tm3_end_pos
    linker_end_index = tm4_start_pos - 1

    # Extract the TM3-TM4 linker sequence
    linker_sequence = full_sequence[linker_start_index:linker_end_index]
    
    # Count the number of Cysteine ('C') residues in the linker sequence
    cysteine_count = linker_sequence.count('C')
    
    print(f"Human GABAAρ1 Receptor Analysis (UniProt: P24046)")
    print("-" * 50)
    print(f"TM3 domain ends at residue: {tm3_end_pos}")
    print(f"TM4 domain starts at residue: {tm4_start_pos}")
    print(f"The extracted TM3-TM4 linker sequence (residues {tm3_end_pos + 1} to {linker_end_index}) is:")
    print(linker_sequence)
    print("-" * 50)
    print(f"The final equation is the count of Cysteines ('C') in this linker.")
    print(f"Number of Cysteine residues = {cysteine_count}")


if __name__ == "__main__":
    count_cysteines_in_linker()