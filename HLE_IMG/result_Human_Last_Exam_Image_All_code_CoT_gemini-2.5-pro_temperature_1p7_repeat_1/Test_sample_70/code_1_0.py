def solve_reaction():
    """
    Identifies the two pericyclic reactions in the given thermal transformation.
    """
    
    reaction_1_electrons = 4
    reaction_1_type = "electrocyclic ring-opening"
    reaction_1_stereochem = "conrotatory"

    reaction_2_electrons = 6
    reaction_2_type = "electrocyclic ring-closure"
    reaction_2_stereochem = "disrotatory"

    print("The thermal transformation involves two sequential pericyclic reactions:")
    print("\nReaction 1:")
    print(f"- Type: A {reaction_1_stereochem} {reaction_1_electrons}π-electron {reaction_1_type}.")
    print("- Details: The four-membered cyclobutene ring opens to form a ten-membered ring intermediate (cyclodecapentaene).")
    
    print("\nReaction 2:")
    print(f"- Type: A {reaction_2_stereochem} {reaction_2_electrons}π-electron {reaction_2_type}.")
    print("- Details: A 1,3,5-hexatriene portion of the intermediate cyclizes to form the final product containing two fused six-membered rings.")

    # The final concise answer identifying both reactions.
    final_answer = (
        f"The two specific reactions are a {reaction_1_electrons}π-electron electrocyclic ring-opening, "
        f"followed by a {reaction_2_electrons}π-electron electrocyclic ring-closure."
    )
    print(f"\nSummary: {final_answer}")

solve_reaction()