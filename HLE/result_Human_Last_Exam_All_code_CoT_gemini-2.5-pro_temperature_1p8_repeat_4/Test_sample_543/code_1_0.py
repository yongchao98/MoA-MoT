def get_iupac_name():
    """
    This function determines and prints the IUPAC name for the product of the reaction between
    1,3-dibromo-2-iodobenzene and excess phenyl magnesium bromide.

    Reaction analysis:
    1. Reactants: 1,3-dibromo-2-iodobenzene and phenyl magnesium bromide (a Grignard reagent).
    2. Reaction Type: Grignard cross-coupling.
    3. Leaving Group Reactivity: The reactivity of aryl halides in this reaction is I > Br.
    4. Product Formation: The phenyl group from the Grignard reagent will preferentially replace
       the most reactive leaving group, the iodine atom at position 2. The less reactive
       bromine atoms at positions 1 and 3 will remain.
    5. Final Product Structure: The benzene ring will be substituted with bromo groups at positions 1 and 3,
       and a phenyl group at position 2.
    """
    # IUPAC Name Construction:
    # Parent: benzene
    # Substituents: 1-bromo, 3-bromo, 2-phenyl
    # Alphabetical order: bromo, phenyl
    # Final Name: 1,3-dibromo-2-phenylbenzene
    
    # Numbers in the final name
    position_1 = 1
    position_3 = 3
    position_2 = 2
    
    # Constructing the name string. The prompt asks to output each number in the final equation.
    # We will interpret this as printing the final name which contains the necessary numbers.
    
    name = f"{position_1},{position_3}-dibromo-{position_2}-phenylbenzene"
    
    print(f"The IUPAC name of the product is: {name}")

get_iupac_name()