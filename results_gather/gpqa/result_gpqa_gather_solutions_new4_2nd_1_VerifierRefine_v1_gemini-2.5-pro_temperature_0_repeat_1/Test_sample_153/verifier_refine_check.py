import re

def check_chemistry_answer():
    """
    This function checks the correctness of the final answer for a chemical identification problem
    based on provided spectral data (MS, IR, 1H NMR).
    """

    # --- Define the spectral data from the question ---
    spectral_data = {
        "ms": {
            "m_z": 156,
            "m_z_plus_2_ratio": 0.32,  # Approx 1/3, indicating one Cl atom
        },
        "ir": {
            "oh_stretch": "broad 3500-2700 cm^-1",
            "co_stretch": 1720,
        },
        "nmr": {
            "acid_h_shift": 11.0,
            "aromatic_protons": 4,
            "aromatic_pattern": "two doublets",
        }
    }

    # --- Define the properties of the candidate molecules ---
    # These properties will be checked against the spectral data.
    candidates = {
        'A': {
            'name': '3-Chloro-2-hydroxybenzaldehyde',
            'formula': 'C7H5ClO2',
            'is_carboxylic_acid': False,
            'substitution': '1,2,3-trisubstituted',
            'aromatic_protons': 3
        },
        'B': {
            'name': '2-chlorobenzoic acid',
            'formula': 'C7H5ClO2',
            'is_carboxylic_acid': True,
            'substitution': 'ortho', # 1,2-disubstituted
            'aromatic_protons': 4
        },
        'C': {
            'name': '4-chlorobenzoic acid',
            'formula': 'C7H5ClO2',
            'is_carboxylic_acid': True,
            'substitution': 'para', # 1,4-disubstituted
            'aromatic_protons': 4
        },
        'D': {
            'name': 'Phenyl chloroformate',
            'formula': 'C7H5ClO2',
            'is_carboxylic_acid': False,
            'substitution': 'monosubstituted',
            'aromatic_protons': 5
        }
    }

    # --- The final answer provided by the LLM ---
    final_answer_key = 'C'
    
    # --- Verification Logic ---
    
    chosen_candidate = candidates.get(final_answer_key)

    if not chosen_candidate:
        return f"Invalid answer key '{final_answer_key}'. The options are A, B, C, D."

    # 1. Check for Carboxylic Acid Functional Group
    # The IR data (broad peak 3500-2700 cm⁻¹) and NMR data (11.0 ppm signal) are
    # definitive evidence for a carboxylic acid.
    if not chosen_candidate['is_carboxylic_acid']:
        return (f"Incorrect: The answer '{chosen_candidate['name']}' is wrong. "
                f"The IR data (broad peak 3500-2700 cm⁻¹) and the ¹H NMR signal at 11.0 ppm "
                f"are definitive evidence for a carboxylic acid functional group, which "
                f"'{chosen_candidate['name']}' does not have.")

    # 2. Check the ¹H NMR Aromatic Region
    # The data shows 4 aromatic protons (2H + 2H).
    if chosen_candidate['aromatic_protons'] != spectral_data['nmr']['aromatic_protons']:
        return (f"Incorrect: The answer '{chosen_candidate['name']}' is wrong. "
                f"The ¹H NMR data shows a total of {spectral_data['nmr']['aromatic_protons']} aromatic protons (2H + 2H). "
                f"However, the proposed structure '{chosen_candidate['name']}' has {chosen_candidate['aromatic_protons']} aromatic protons.")

    # 3. Check the ¹H NMR Aromatic Splitting Pattern
    # The data shows two doublets, which is a classic pattern for a para-substituted ring.
    if chosen_candidate['substitution'] != 'para':
        return (f"Incorrect: The answer '{chosen_candidate['name']}' is wrong. "
                f"The ¹H NMR aromatic region shows two doublets, each for 2 protons. This symmetrical pattern "
                f"is characteristic of a para (1,4) disubstituted benzene ring. "
                f"The proposed structure has '{chosen_candidate['substitution']}' substitution, which would produce a different, more complex NMR spectrum.")

    # 4. Check Mass Spectrometry (as a final confirmation)
    # This check is somewhat redundant as all options have the same formula, but it's good practice.
    formula = chosen_candidate['formula']
    atom_mass = {'C': 12, 'H': 1, 'O': 16, 'Cl': 35} # Nominal mass with most abundant isotopes
    mass = sum(atom_mass.get(atom, 0) * int(count or 1) for atom, count in re.findall(r'([A-Z][a-z]*)(\d*)', formula))
    
    if mass != spectral_data['ms']['m_z']:
        return (f"Incorrect: The answer '{chosen_candidate['name']}' is wrong. "
                f"The mass spectrum shows a molecular ion peak at m/z = {spectral_data['ms']['m_z']}. "
                f"The nominal mass of the proposed structure is {mass}.")
    
    if 'Cl' not in formula:
        return (f"Incorrect: The answer '{chosen_candidate['name']}' is wrong. "
                f"The mass spectrum's M+2 peak at 32% intensity indicates the presence of a chlorine atom, "
                f"which is missing from the proposed structure's formula.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_chemistry_answer()
print(result)