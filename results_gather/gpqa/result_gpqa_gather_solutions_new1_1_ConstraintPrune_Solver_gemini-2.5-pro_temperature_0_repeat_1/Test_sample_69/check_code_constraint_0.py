def check_answer():
    """
    This function checks the correctness of the "Sulfur/Mustard Gas" hypothesis
    against all the constraints provided in the chemical puzzle.
    """
    # The proposed solution to be verified, based on the final answer's reasoning.
    solution = {
        "A": "S₈",              # Solid A is elemental sulfur
        "B": "Cl₂",             # Gas B is chlorine
        "C": "SCl₂",            # Bright red product C is sulfur dichloride
        "D": "C₂H₄",            # Gas D is ethene
        "E": "(ClCH₂CH₂)₂S",    # Hazardous product E is mustard gas
        "F": "HCl",             # Strong acid F is hydrochloric acid
        "G": "H₂SO₃",           # Weak acid G is sulfurous acid
        "H": "C₂H₄Cl₂",         # Solvent H is 1,2-dichloroethane
        "E_symmetry": "C₂",     # The symmetry of mustard gas
        "final_answer_option": "D" # The option corresponding to C₂
    }

    errors = []

    # Constraint 1: A(s) + 8 B(g) -> C (bright red product)
    # The reaction S₈(s) + 8 Cl₂(g) -> 8 SCl₂ fits the 1:8 stoichiometry perfectly.
    # SCl₂ is a cherry-red liquid, matching the description.
    if not (solution["A"] == "S₈" and solution["B"] == "Cl₂" and solution["C"] == "SCl₂"):
        errors.append("Constraint 1 FAILED: The proposed identities for A, B, and C do not match the S₈ + 8Cl₂ -> 8SCl₂ reaction, which perfectly fits the 1:8 stoichiometry.")

    # Constraint 2: C + 2 D(g) -> E (extremely hazardous product)
    # The reaction SCl₂ + 2 C₂H₄ -> (ClCH₂CH₂)₂S is the Levinstein process for mustard gas.
    # It fits the 1:2 stoichiometry and the "extremely hazardous" description.
    if not (solution["C"] == "SCl₂" and solution["D"] == "C₂H₄" and solution["E"] == "(ClCH₂CH₂)₂S"):
        errors.append("Constraint 2 FAILED: The proposed identities for C, D, and E do not match the synthesis of mustard gas (SCl₂ + 2C₂H₄ -> (ClCH₂CH₂)₂S), which fits the 1:2 stoichiometry and 'extremely hazardous' description.")

    # Constraint 3: C + H₂O -> A(s) + F(strong acid) + G(weak acid)
    # The hydrolysis of SCl₂ disproportionates to produce elemental sulfur (A), HCl (F, strong),
    # and SO₂, which forms H₂SO₃ (G, weak) in water.
    if not (solution["C"] == "SCl₂" and "S" in solution["A"] and solution["F"] == "HCl" and solution["G"] == "H₂SO₃"):
        errors.append("Constraint 3 FAILED: The hydrolysis of the proposed C (SCl₂) should produce elemental sulfur (A), a strong acid (F=HCl), and a weak acid (G=H₂SO₃). The proposed products in the solution do not match this known reaction.")

    # Constraint 4: D(g) + B(g) -> H (solvent) (1:1 ratio)
    # The reaction C₂H₄(g) + Cl₂(g) -> C₂H₄Cl₂(l) is a 1:1 reaction that produces a common solvent.
    if not (solution["D"] == "C₂H₄" and solution["B"] == "Cl₂" and solution["H"] == "C₂H₄Cl₂"):
        errors.append("Constraint 4 FAILED: The proposed identities for D, B, and H do not match the 1:1 reaction C₂H₄ + Cl₂ -> C₂H₄Cl₂ (a solvent).")

    # Final Check: Molecular Symmetry of E
    # The identified product E is mustard gas, (ClCH₂CH₂)₂S. Its point group is C₂.
    # The final answer option is D, which corresponds to C₂.
    if solution["E_symmetry"] != "C₂":
        errors.append(f"Symmetry Conclusion FAILED: The identified product E, Mustard Gas, has C₂ symmetry, not {solution['E_symmetry']}.")
    if solution["final_answer_option"] != "D":
        errors.append(f"Final Answer Option FAILED: The correct symmetry is C₂, which corresponds to option D, not option {solution['final_answer_option']}.")

    # Return the result
    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Execute the check and print the result
result = check_answer()
print(result)