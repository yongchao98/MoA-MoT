import math

def check_answer():
    """
    This function checks the correctness of the provided answer by codifying the spectral data
    and the properties of the candidate molecules.
    """

    # 1. Define the experimental data from the question
    experimental_data = {
        "mass_spec": {
            "m_z": 156,
            "m_z_plus_2_ratio": 0.32  # approx 1/3, indicating one Cl atom
        },
        "ir": {
            "has_broad_oh_3500_2700": True,  # Characteristic of a carboxylic acid
            "has_co_1720": True             # Characteristic of a conjugated carbonyl
        },
        "nmr": {
            "has_cooh_proton": True,        # Signal at 11.0 ppm
            "aromatic_protons": 4,
            "aromatic_pattern": "para"      # Two doublets, 2H each, indicates 1,4-disubstitution
        }
    }

    # 2. Define the properties of the candidate molecules
    candidates = {
        "A": {
            "name": "4-chlorobenzoic acid",
            "mw": 156, # C7H5(35)ClO2 = 7*12 + 5*1 + 35 + 2*16 = 156
            "has_chlorine": True,
            "functional_group": "carboxylic acid",
            "nmr_aromatic_protons": 4,
            "nmr_aromatic_pattern": "para"
        },
        "B": {
            "name": "2-chlorobenzoic acid",
            "mw": 156,
            "has_chlorine": True,
            "functional_group": "carboxylic acid",
            "nmr_aromatic_protons": 4,
            "nmr_aromatic_pattern": "ortho" # Asymmetrical, would give complex multiplets
        },
        "C": {
            "name": "Phenyl chloroformate",
            "mw": 156,
            "has_chlorine": True,
            "functional_group": "acid chloride/ester", # Not a carboxylic acid
            "nmr_aromatic_protons": 5, # Monosubstituted ring
            "nmr_aromatic_pattern": "mono"
        },
        "D": {
            "name": "3-Chloro-2-hydroxybenzaldehyde",
            "mw": 156,
            "has_chlorine": True,
            "functional_group": "aldehyde/phenol", # Not a carboxylic acid
            "nmr_aromatic_protons": 3, # Trisubstituted ring
            "nmr_aromatic_pattern": "tri"
        }
    }

    llm_answer_choice = "A"
    llm_answer_properties = candidates[llm_answer_choice]

    # 3. Perform the checks
    # Check Mass Spec
    if llm_answer_properties["mw"] != experimental_data["mass_spec"]["m_z"]:
        return f"Incorrect. The molecular weight of {llm_answer_properties['name']} ({llm_answer_properties['mw']}) does not match the experimental M+ peak at m/z = {experimental_data['mass_spec']['m_z']}."
    if not llm_answer_properties["has_chlorine"]:
        return f"Incorrect. The mass spectrum's M+2 peak indicates a chlorine atom is present, but {llm_answer_properties['name']} does not contain chlorine."

    # Check IR
    if llm_answer_properties["functional_group"] != "carboxylic acid" and experimental_data["ir"]["has_broad_oh_3500_2700"]:
        return f"Incorrect. The IR spectrum strongly indicates a carboxylic acid, but {llm_answer_properties['name']} is a(n) {llm_answer_properties['functional_group']}."

    # Check NMR
    if llm_answer_properties["functional_group"] != "carboxylic acid" and experimental_data["nmr"]["has_cooh_proton"]:
        return f"Incorrect. The NMR spectrum shows a carboxylic acid proton at 11.0 ppm, which is absent in {llm_answer_properties['name']}."
    if llm_answer_properties["nmr_aromatic_protons"] != experimental_data["nmr"]["aromatic_protons"]:
        return f"Incorrect. The NMR shows {experimental_data['nmr']['aromatic_protons']} aromatic protons, but {llm_answer_properties['name']} has {llm_answer_properties['nmr_aromatic_protons']}."
    if llm_answer_properties["nmr_aromatic_pattern"] != experimental_data["nmr"]["aromatic_pattern"]:
        return f"Incorrect. The NMR splitting pattern (two doublets) indicates a symmetrical '{experimental_data['nmr']['aromatic_pattern']}' substitution, but {llm_answer_properties['name']} has a '{llm_answer_properties['nmr_aromatic_pattern']}' substitution pattern, which would produce a different signal."

    # If all checks pass for the given answer, we verify that no other answer could be correct.
    for choice, props in candidates.items():
        if choice == llm_answer_choice:
            continue
        # Check if any other candidate could also be correct
        if (props["mw"] == experimental_data["mass_spec"]["m_z"] and
            props["functional_group"] == "carboxylic acid" and
            props["nmr_aromatic_protons"] == experimental_data["nmr"]["aromatic_protons"] and
            props["nmr_aromatic_pattern"] == experimental_data["nmr"]["aromatic_pattern"]):
            return f"Incorrect. The answer {llm_answer_choice} is plausible, but option {choice} ({props['name']}) also fits all criteria. The analysis is ambiguous."

    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)