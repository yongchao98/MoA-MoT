def check_spectroscopic_data(candidate_letter):
    """
    Checks a candidate structure against the provided spectroscopic data.

    Args:
        candidate_letter (str): The letter of the candidate to check ('A', 'B', 'C', or 'D').

    Returns:
        str: "Correct" if the candidate matches all data, otherwise a string
             explaining the reason for the mismatch.
    """
    # --- Define Candidate Properties ---
    candidates = {
        "A": {
            "name": "ethyl 4-aminobenzoate",
            "amine_type": "primary",  # Has -NH2 group
            "carbonyl_type": "ester", # Has -COO- group
            "aromatic_pattern": "para", # 1,4-disubstituted
            "alkyl_group_nmr_key_shift": "quartet_at_4.5" # -O-CH2CH3
        },
        "B": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "amine_type": "secondary_amide", # Has -NH- group
            "carbonyl_type": "amide",
            "aromatic_pattern": "para",
            "alkyl_group_nmr_key_shift": "quartet_at_4.0" # -O-CH2CH3
        },
        "C": {
            "name": "4-aminophenyl propionate",
            "amine_type": "primary",
            "carbonyl_type": "ester",
            "aromatic_pattern": "para",
            "alkyl_group_nmr_key_shift": "quartet_at_2.5" # -CO-CH2CH3
        },
        "D": {
            "name": "3-ethoxybenzamide",
            "amine_type": "primary_amide", # Has -CONH2 group
            "carbonyl_type": "amide",
            "aromatic_pattern": "meta", # 1,3-disubstituted
            "alkyl_group_nmr_key_shift": "quartet_at_4.0" # -O-CH2CH3
        }
    }

    if candidate_letter not in candidates:
        return f"Invalid candidate letter '{candidate_letter}'. Please choose from 'A', 'B', 'C', or 'D'."

    candidate = candidates[candidate_letter]
    name = candidate["name"]

    # --- Constraint 1: IR Spectrum (Amine Type) ---
    # Data: Two N-H bands at 3420/3325 cm-1 indicate a primary amine (-NH2).
    if candidate["amine_type"] not in ["primary", "primary_amide"]:
        return (f"Incorrect. Candidate {candidate_letter} ({name}) is a {candidate['amine_type'].replace('_', ' ')}, "
                f"but the IR spectrum (two bands at 3420 and 3325 cm-1) is characteristic of a primary amine.")

    # --- Constraint 2: IR Spectrum (Carbonyl Type) ---
    # Data: Strong C=O band at 1720 cm-1 indicates a conjugated ester.
    if candidate["carbonyl_type"] == "amide":
        return (f"Incorrect. Candidate {candidate_letter} ({name}) is an amide. The IR C=O stretch at 1720 cm-1 "
                f"is more consistent with a conjugated ester than an amide (typically < 1700 cm-1).")

    # --- Constraint 3: 1H NMR (Aromatic Pattern) ---
    # Data: Two doublets at 7.0 and 8.0 ppm indicate a 1,4- (para) substituted ring.
    if candidate["aromatic_pattern"] != "para":
        return (f"Incorrect. Candidate {candidate_letter} ({name}) is {candidate['aromatic_pattern']}-substituted, "
                f"but the NMR spectrum (two doublets at 7.0 and 8.0 ppm) indicates a para-substitution pattern.")

    # --- Constraint 4: 1H NMR (Alkyl Group) ---
    # Data: Quartet at 4.5 ppm indicates an ethyl group attached to an oxygen (-O-CH2CH3).
    if candidate["alkyl_group_nmr_key_shift"] != "quartet_at_4.5":
        if candidate["alkyl_group_nmr_key_shift"] == "quartet_at_2.5":
             return (f"Incorrect. Candidate {candidate_letter} ({name}) has a propionate group (-CO-CH2CH3). "
                     f"The NMR quartet for this -CH2- group would be around 2.5 ppm, not 4.5 ppm as observed.")
        else:
             return (f"Incorrect. The NMR data for candidate {candidate_letter} ({name}) is inconsistent with the observed spectrum.")


    # If all constraints are satisfied, the answer is correct.
    return "Correct"

# Example usage of the checking code:
# print(f"Checking A: {check_spectroscopic_data('A')}")
# print(f"Checking B: {check_spectroscopic_data('B')}")
# print(f"Checking C: {check_spectroscopic_data('C')}")
# print(f"Checking D: {check_spectroscopic_data('D')}")