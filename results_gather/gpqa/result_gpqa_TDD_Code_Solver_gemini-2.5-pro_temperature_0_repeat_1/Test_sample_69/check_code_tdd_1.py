def check_answer():
    """
    Checks the correctness of the LLM's answer by verifying its chemical reasoning.
    
    The code checks the following logical steps from the reasoning:
    1. The identity of the hazardous product E is [ClF4]-.
    2. The molecular geometry of [ClF4]- is square planar.
    3. The point group of a square planar molecule is D4h.
    4. The final answer 'D' corresponds to D4h.
    """

    # A knowledge base of molecular geometries and point groups
    # This represents established chemical principles (like VSEPR and Group Theory)
    molecular_data = {
        '[ClF4]-': {
            'vsepr_type': 'AX4E2',
            'shape': 'square planar',
            'point_group': 'D4h'
        },
        'H2O': {
            'vsepr_type': 'AX2E2',
            'shape': 'bent',
            'point_group': 'C2v'
        },
        'CO2': {
            'vsepr_type': 'AX2',
            'shape': 'linear',
            'point_group': 'Dinfh' # Represents D∞h
        },
        'CH2Cl2': {
            'vsepr_type': 'AX4',
            'shape': 'tetrahedral',
            'point_group': 'C2v'
        }
    }
    
    # Mapping of options to point groups from the question
    answer_options = {
        'A': 'C2v',
        'B': 'Dinfh',
        'C': 'C2',
        'D': 'D4h'
    }

    # --- Verification based on the LLM's reasoning ---

    # 1. LLM identifies the hazardous product E as [ClF4]-.
    identified_molecule_E = '[ClF4]-'

    # 2. LLM claims E has a square planar geometry. Let's check our data.
    if identified_molecule_E not in molecular_data:
        return f"Reasoning Error: The molecule {identified_molecule_E} identified by the LLM is not in our verification database."

    expected_shape = 'square planar'
    actual_shape = molecular_data[identified_molecule_E]['shape']
    if actual_shape != expected_shape:
        return f"Reasoning Error: The LLM claims the shape of {identified_molecule_E} is {expected_shape}, but it is actually {actual_shape}."

    # 3. LLM claims a square planar molecule has D4h symmetry. Let's check our data.
    expected_point_group = 'D4h'
    actual_point_group = molecular_data[identified_molecule_E]['point_group']
    if actual_point_group != expected_point_group:
        return f"Reasoning Error: The LLM claims the point group for a {actual_shape} molecule like {identified_molecule_E} is {expected_point_group}, but it is {actual_point_group}."

    # 4. LLM concludes the answer is D. Let's check if D matches the derived point group.
    llm_final_answer_option = 'D'
    if answer_options.get(llm_final_answer_option) != actual_point_group:
        return f"Answer Incorrect: The reasoning correctly leads to the point group {actual_point_group}, but the final answer provided is {llm_final_answer_option} ({answer_options.get(llm_final_answer_option)}), which is a mismatch."

    # If all checks on the reasoning pass and lead to the given answer, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)