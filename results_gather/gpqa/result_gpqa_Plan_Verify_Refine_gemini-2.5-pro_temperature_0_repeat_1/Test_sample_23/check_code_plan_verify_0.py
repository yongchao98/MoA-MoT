import math

def check_correctness_of_llm_answer():
    """
    This function checks the correctness of the LLM's reasoning and final answer
    by verifying that the proposed chemical compounds satisfy all constraints
    given in the problem description.
    """

    # --- Data Definitions based on LLM's reasoning ---

    # Using precise atomic masses for better accuracy
    ATOMIC_MASS = {'C': 12.011, 'H': 1.008}

    # Define the properties of the proposed compounds
    # 'decolorizes_br2_water' is a proxy for being unsaturated (non-aromatic)
    # 'carbon_skeleton_id' checks that hydrogenation leads to the same base molecule
    compounds = {
        "cyclohexane": {
            "formula": (6, 12), "decolorizes_br2_water": False,
            "is_saturated": True, "has_conjugated_bonds": False,
            "carbon_skeleton_id": "C6_cyclic"
        },
        "benzene": {
            "formula": (6, 6), "decolorizes_br2_water": False,
            "is_saturated": False, "has_conjugated_bonds": True,
            "carbon_skeleton_id": "C6_cyclic"
        },
        "cyclohexene": {
            "formula": (6, 10), "decolorizes_br2_water": True,
            "is_saturated": False, "has_conjugated_bonds": False,
            "carbon_skeleton_id": "C6_cyclic"
        },
        "1,4-cyclohexadiene": {
            "formula": (6, 8), "decolorizes_br2_water": True,
            "is_saturated": False, "has_conjugated_bonds": False,
            "carbon_skeleton_id": "C6_cyclic"
        }
    }

    # The LLM's proposed solution identifies the compounds as:
    Z_name = "cyclohexane"
    Y_components_names = ["benzene", "cyclohexane"]
    X_components_names = ["1,4-cyclohexadiene", "cyclohexene"]
    llm_final_answer = 18  # The total number of H atoms in X from option A

    # --- Verification Steps ---

    error_messages = []

    # 1. Verify Substance Z (Cyclohexane)
    Z = compounds[Z_name]
    c_atoms, h_atoms = Z["formula"]
    h_mass_fraction = (h_atoms * ATOMIC_MASS['H']) / (c_atoms * ATOMIC_MASS['C'] + h_atoms * ATOMIC_MASS['H'])
    if not math.isclose(h_mass_fraction, 0.1428, rel_tol=0.015):
        error_messages.append(f"Constraint Fail (Z): Calculated H mass fraction for {Z_name} is {h_mass_fraction:.4f}, not ~14.28%.")
    if not Z["is_saturated"]:
        error_messages.append(f"Constraint Fail (Z): Proposed substance {Z_name} is not saturated.")

    # 2. Verify Mixture Y
    Y1, Y2 = compounds[Y_components_names[0]], compounds[Y_components_names[1]]
    if Y1["decolorizes_br2_water"] or Y2["decolorizes_br2_water"]:
        error_messages.append("Constraint Fail (Y): At least one component of proposed mixture Y decolorizes bromine water.")
    if Y1["carbon_skeleton_id"] != Z["carbon_skeleton_id"] or Y2["carbon_skeleton_id"] != Z["carbon_skeleton_id"]:
        error_messages.append("Constraint Fail (Y): Hydrogenation of proposed mixture Y does not yield only Z.")

    # 3. Verify Mixture X
    X1, X2 = compounds[X_components_names[0]], compounds[X_components_names[1]]
    if not (X1["decolorizes_br2_water"] and X2["decolorizes_br2_water"]):
        error_messages.append("Constraint Fail (X): At least one component of proposed mixture X does not decolorize bromine water.")
    if X1["has_conjugated_bonds"] or X2["has_conjugated_bonds"]:
        error_messages.append("Constraint Fail (X): At least one component of proposed mixture X has conjugated bonds.")
    if X1["carbon_skeleton_id"] != Z["carbon_skeleton_id"] or X2["carbon_skeleton_id"] != Z["carbon_skeleton_id"]:
        error_messages.append("Constraint Fail (X): Hydrogenation of proposed mixture X does not yield only Z.")

    # 4. Verify Disproportionation Reaction (X -> Y)
    x_total_h = X1["formula"][1] + X2["formula"][1]
    y_total_h = Y1["formula"][1] + Y2["formula"][1]
    x_total_c = X1["formula"][0] + X2["formula"][0]
    y_total_c = Y1["formula"][0] + Y2["formula"][0]
    if x_total_c != y_total_c or x_total_h != y_total_h:
        error_messages.append(f"Constraint Fail (Reaction): The disproportionation X -> Y is not balanced. X has C{x_total_c}H{x_total_h}, Y has C{y_total_c}H{y_total_h}.")

    # 5. Verify Final Answer Calculation
    calculated_h_atoms_in_X = X1["formula"][1] + X2["formula"][1]
    if calculated_h_atoms_in_X != llm_final_answer:
        error_messages.append(f"Final Answer Mismatch: The total H atoms in the proposed mixture X is {calculated_h_atoms_in_X}, but the answer is {llm_final_answer}.")

    # --- Return Result ---
    if not error_messages:
        return "Correct"
    else:
        return "\n".join(error_messages)

# Execute the check
result = check_correctness_of_llm_answer()
print(result)