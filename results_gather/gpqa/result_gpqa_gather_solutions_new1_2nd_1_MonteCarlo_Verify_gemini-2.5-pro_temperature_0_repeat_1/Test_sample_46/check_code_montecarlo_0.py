import re

def check_spectroscopic_data():
    """
    Checks if the proposed compound (ethyl 4-aminobenzoate) matches the given spectral data.
    """
    # --- Data from the Question ---
    spectral_data = {
        "ir_nh_bands": 2,  # Two bands at 3420 and 3325 cm-1 imply a primary amine.
        "ir_co_freq": 1720, # A frequency typical for a conjugated ester.
        "nmr_aromatic_pattern": "para", # Two doublets, each 2H, implies 1,4-substitution.
        "nmr_quartet_ppm": 4.5, # The chemical shift of the -CH2- quartet.
        "nmr_amine_protons": True # A signal at 4.0 ppm for 2H confirms a primary amine.
    }

    # --- Candidate Compounds and their Expected Properties ---
    # This dictionary serves as our chemical knowledge base.
    candidates = {
        "A": {
            "name": "4-aminophenyl propionate",
            "amine_type": "primary",
            "carbonyl_type": "ester",
            "substitution": "para",
            "quartet_env": "propionyl", # -CO-CH2-
            "expected_quartet_ppm_range": (2.2, 2.8)
        },
        "B": {
            "name": "ethyl 4-aminobenzoate",
            "amine_type": "primary",
            "carbonyl_type": "ester",
            "substitution": "para",
            "quartet_env": "ethyl_ester", # -O-CH2-
            "expected_quartet_ppm_range": (4.1, 4.6)
        },
        "C": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "amine_type": "secondary_amide", # Not a primary amine
            "carbonyl_type": "amide",
            "substitution": "para",
            "quartet_env": "ethoxy", # Ar-O-CH2-
            "expected_quartet_ppm_range": (3.9, 4.2)
        },
        "D": {
            "name": "3-ethoxybenzamide",
            "amine_type": "primary_amide",
            "carbonyl_type": "amide",
            "substitution": "meta", # Not para
            "quartet_env": "ethoxy", # Ar-O-CH2-
            "expected_quartet_ppm_range": (3.9, 4.2)
        }
    }

    # The final answer to check
    final_answer_key = "B"
    chosen_candidate = candidates[final_answer_key]

    # --- Verification Steps ---

    # 1. Check Amine Type (from IR N-H bands)
    if spectral_data["ir_nh_bands"] == 2:
        if chosen_candidate["amine_type"] not in ["primary", "primary_amide"]:
            return f"Incorrect. The IR data shows two N-H bands, indicating a primary amine (-NH2). However, {chosen_candidate['name']} is a {chosen_candidate['amine_type']}."

    # 2. Check Carbonyl Type (from IR C=O frequency)
    # Conjugated esters are ~1715-1730 cm-1. Conjugated amides are typically lower, ~1650-1690 cm-1.
    if spectral_data["ir_co_freq"] > 1700 and chosen_candidate["carbonyl_type"] != "ester":
        return f"Incorrect. The IR C=O frequency of {spectral_data['ir_co_freq']} cm-1 is characteristic of an ester, but {chosen_candidate['name']} is an {chosen_candidate['carbonyl_type']}."

    # 3. Check Aromatic Substitution Pattern (from NMR)
    if spectral_data["nmr_aromatic_pattern"] == "para" and chosen_candidate["substitution"] != "para":
        return f"Incorrect. The NMR data shows a para-substitution pattern, but {chosen_candidate['name']} is {chosen_candidate['substitution']}-substituted."

    # 4. Check the Decisive NMR Quartet Chemical Shift
    quartet_ppm = spectral_data["nmr_quartet_ppm"]
    expected_range = chosen_candidate["expected_quartet_ppm_range"]
    if not (expected_range[0] <= quartet_ppm <= expected_range[1]):
        return (f"Incorrect. The most decisive piece of evidence is the NMR quartet at {quartet_ppm} ppm. "
                f"This indicates an ethyl group attached to an ester oxygen (-O-CH2-). "
                f"For {chosen_candidate['name']}, the quartet from its '{chosen_candidate['quartet_env']}' group is expected in the range {expected_range[0]}-{expected_range[1]} ppm, which does not match the data.")

    # If all checks pass for the chosen candidate, the answer is correct.
    return "Correct"

# Run the verification
result = check_spectroscopic_data()
print(result)