import collections

def check_spectroscopy_answer():
    """
    Checks the correctness of the identified compound based on spectral data.
    """
    question_data = {
        "formula": "C9H11NO2",
        "IR": {
            "N-H_stretches": 2,  # Two bands at 3420, 3325 cm-1 indicate a primary amine
            "C=O_stretch": 1720  # cm-1
        },
        "NMR": {
            "signals": [
                {"ppm": 1.20, "type": "t", "H": 3},
                {"ppm": 4.0, "type": "bs", "H": 2},
                {"ppm": 4.5, "type": "q", "H": 2},
                {"ppm": 7.0, "type": "d", "H": 2},
                {"ppm": 8.0, "type": "d", "H": 2}
            ]
        }
    }

    # Properties of the candidate molecules
    candidates = {
        "A": {
            "name": "3-ethoxybenzamide",
            "formula": "C9H11NO2",
            "functional_groups": ["primary_amide", "ether", "meta_disubstituted_ring"],
            "ir_c_double_o_range": (1650, 1690),
            "nmr_aromatic_pattern": "meta" # Expects 3-4 signals, not two doublets
        },
        "B": {
            "name": "ethyl 4-aminobenzoate",
            "formula": "C9H11NO2",
            "functional_groups": ["primary_amine", "ester", "para_disubstituted_ring"],
            "ir_c_double_o_range": (1715, 1730), # Conjugated ester
            "nmr_aromatic_pattern": "para", # Expects two doublets
            "has_ethyl_ester": True
        },
        "C": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "formula": "C9H11NO2",
            "functional_groups": ["secondary_amide", "ether", "para_disubstituted_ring"],
            "ir_c_double_o_range": (1650, 1680),
            "nmr_aromatic_pattern": "para"
        },
        "D": {
            "name": "4-aminophenyl propionate",
            "formula": "C9H11NO2",
            "functional_groups": ["primary_amine", "ester", "para_disubstituted_ring"],
            "ir_c_double_o_range": (1735, 1755), # Phenyl ester
            "nmr_aromatic_pattern": "para",
            "has_ethyl_ester": False # Has a propionate group, not an ethyl ester
        }
    }

    llm_answer = "B"
    
    # --- Verification Logic ---
    
    # Check the proposed correct answer first
    correct_candidate = candidates[llm_answer]
    
    # 1. Check Molecular Formula
    if correct_candidate["formula"] != question_data["formula"]:
        return f"Incorrect. The formula for {correct_candidate['name']} is not {question_data['formula']}."

    # 2. Check IR data
    # Check for primary amine
    if "primary_amine" not in correct_candidate["functional_groups"]:
        return f"Incorrect. The IR shows two N-H bands (3420, 3325 cm-1), characteristic of a primary amine (-NH2), which {correct_candidate['name']} does not have."
    # Check C=O stretch
    c_o_range = correct_candidate["ir_c_double_o_range"]
    if not (c_o_range[0] <= question_data["IR"]["C=O_stretch"] <= c_o_range[1]):
        return f"Incorrect. The IR C=O stretch at 1720 cm-1 is not in the expected range for the carbonyl in {correct_candidate['name']} ({c_o_range[0]}-{c_o_range[1]} cm-1)."

    # 3. Check NMR data
    # Check for ethyl group attached to oxygen (explains t @ 1.2 and q @ 4.5)
    if not correct_candidate.get("has_ethyl_ester", False):
         return f"Incorrect. The NMR shows a triplet at 1.2 ppm and a quartet at 4.5 ppm, indicating an ethyl group attached to an oxygen (-OCH2CH3). This is not present in {correct_candidate['name']}."
    
    # Check for para-disubstituted ring (explains two doublets in aromatic region)
    if correct_candidate["nmr_aromatic_pattern"] != "para":
        return f"Incorrect. The NMR shows two doublets in the aromatic region, indicating a 1,4-(para) substituted ring. {correct_candidate['name']} has a {correct_candidate['nmr_aromatic_pattern']} substitution pattern."

    # Check if other candidates can be ruled out
    for option, candidate_data in candidates.items():
        if option == llm_answer:
            continue

        # Rule out A: meta substitution and amide C=O
        if option == "A":
            if candidate_data["nmr_aromatic_pattern"] != "para":
                pass # Correctly ruled out
            else:
                return "Logic Error: Candidate A should be meta."
            if not (candidate_data["ir_c_double_o_range"][0] <= 1720 <= candidate_data["ir_c_double_o_range"][1]):
                 pass # Correctly ruled out
            else:
                return "Logic Error: Candidate A's C=O range is wrong."

        # Rule out C: secondary amide (only 1 N-H stretch)
        if option == "C":
            if "secondary_amide" in candidate_data["functional_groups"]:
                pass # Correctly ruled out
            else:
                return "Logic Error: Candidate C should be a secondary amide."

        # Rule out D: propionate group (NMR ethyl signals would be different)
        if option == "D":
            if not candidate_data.get("has_ethyl_ester", True):
                 pass # Correctly ruled out
            else:
                return "Logic Error: Candidate D has a propionate, not an ethyl ester."

    # If all checks pass for the correct answer and others are ruled out, the answer is correct.
    return "Correct"

# Run the check
result = check_spectroscopy_answer()
print(result)