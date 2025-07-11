def count_cysteines_in_linker():
    """
    This script calculates the number of Cysteine residues in the TM3-TM4 linker
    of the human GABAA_rho1 receptor (UniProt P24046).

    The boundaries of the linker are determined from the experimental structure (PDB: 6HUP).
    """

    # Full amino acid sequence for human GABAA_rho1 from UniProt P24046
    full_sequence = "MGFHAGSKLLGFLPLLGAALALGSVPAGSTRYIQQSAPDTLDINLMLDGSVVELTESFHERETVNDYLISRIYLDDRLLLGYDNRLRPGFGGAVTIGTNYTITMTYKCTIDIWLPETFFVNSKKSFIDHDMEYAYEMPCFVLEDGPVLGVTVRLSPSSTGEYIVMTTHFHLKRKIGYFVIQTYLPCIMTVILSQVSFWLNRESVPARTVFGVTTVLTMTTLSISARNSLPKVAYATAMDWFIAVCFAFVFSALIEFATVNYFTKRGYAWDGKSVVPEKPKKVKDPLIKKNQTYIPNLRGPRVLGTPDSRLLKRRNILSLRIKIDRLHIASTVINLPHHPEDGAKMPKKTFETPKNSLEMLKKTSSKKADKNTTMSYIRSKSMFKPDSLTAAIHRNATPSPSQMALEDAAFPRSATATDTYNSIMCGSSVEVPAGLGPSPYPKSPAAKSTASSQASLSSYATAMSRSETVVVSESTGSLPSDLSVMGCCCRASGSYGGGDMLHCTQTEAGEITIETDSFSIVGSKPETPETPA"

    # Define the linker region based on PDB structure 6HUP (1-based indexing)
    # TM3 ends at residue 318
    # TM4 starts at residue 427
    linker_start_residue = 319
    linker_end_residue = 426

    # Convert to 0-based Python slicing indices
    linker_start_index = linker_start_residue - 1
    linker_end_index = linker_end_residue

    # Extract the linker sequence
    linker_sequence = full_sequence[linker_start_index:linker_end_index]

    cysteine_positions = []
    # Find positions of all Cysteine ('C') residues
    for i, amino_acid in enumerate(linker_sequence):
        if amino_acid == 'C':
            # Convert local index back to 1-based residue number in the full sequence
            position = linker_start_index + i + 1
            cysteine_positions.append(position)

    count = len(cysteine_positions)

    print(f"Human GABAA p1 Receptor (UniProt P24046)")
    print(f"TM3-TM4 linker domain is defined as residues {linker_start_residue}-{linker_end_residue} (based on PDB: 6HUP).")
    print("-" * 50)
    
    if count > 0:
        print(f"Found {count} Cysteine residue(s) at position(s): {cysteine_positions}")
        # Build the final "equation" as requested
        equation_parts = ['1'] * count
        equation_str = " + ".join(equation_parts)
        if count > 1:
            equation_str += f" = {count}"
        else:
             equation_str += f" = {count}" # To ensure output for single Cys
        print(f"The final number of Cysteine residues is: {equation_str}")
    else:
        print("No Cysteine residues were found in the TM3-TM4 linker.")
        print("The final number of Cysteine residues is: 0")

count_cysteines_in_linker()
<<<1>>>