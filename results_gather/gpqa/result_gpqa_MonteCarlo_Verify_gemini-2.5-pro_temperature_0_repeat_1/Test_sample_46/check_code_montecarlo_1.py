import collections

def check_spectroscopy_answer():
    """
    This function checks the correctness of the proposed answer for the given spectroscopy problem.
    It analyzes the molecular formula, IR, and 1H NMR data against the properties of each candidate molecule.
    """
    
    # --- Given Data ---
    question_formula = "C9H11NO2"
    ir_data = {
        "N-H_stretches": 2,  # Two bands at 3420 and 3325 cm-1 imply a primary amine/amide
        "C=O_stretch": 1720  # cm-1, characteristic of an ester
    }
    nmr_data = {
        "signals": 5,
        "ethyl_group": {"CH3": {"ppm": 1.20, "type": "t", "H": 3}, "CH2": {"ppm": 4.5, "type": "q", "H": 2}},
        "amine_group": {"ppm": 4.0, "type": "bs", "H": 2},
        "aromatic_pattern": "para-disubstituted" # Two doublets for 2H each
    }
    
    # --- Candidate Structures ---
    candidates = {
        "A": {
            "name": "ethyl 4-aminobenzoate",
            "formula": "C9H11NO2",
            "has_primary_amine": True,
            "carbonyl_type": "ester",  # Expected C=O: ~1715-1735 cm-1
            "has_ethoxy_group": True, # -O-CH2-CH3
            "aromatic_substitution": "para"
        },
        "B": {
            "name": "4-aminophenyl propionate",
            "formula": "C9H11NO2",
            "has_primary_amine": True,
            "carbonyl_type": "ester",
            "has_ethoxy_group": False, # Has -CO-CH2-CH3, not -O-CH2-CH3
            "aromatic_substitution": "para"
        },
        "C": {
            "name": "3-ethoxybenzamide",
            "formula": "C9H11NO2",
            "has_primary_amine": True, # Primary amide has -NH2
            "carbonyl_type": "amide", # Expected C=O: ~1650-1690 cm-1
            "has_ethoxy_group": True,
            "aromatic_substitution": "meta"
        },
        "D": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "formula": "C9H11NO2",
            "has_primary_amine": False, # Has a secondary amide (-NH-)
            "carbonyl_type": "amide",
            "has_ethoxy_group": True,
            "aromatic_substitution": "para"
        }
    }
    
    # The answer provided by the other LLM
    llm_answer = "A"
    
    # --- Verification Logic ---
    
    # Select the candidate corresponding to the given answer
    candidate = candidates[llm_answer]
    
    # 1. Check Molecular Formula
    if candidate["formula"] != question_formula:
        return f"Incorrect. The molecular formula of {candidate['name']} is {candidate['formula']}, which does not match the required {question_formula}."

    # 2. Check IR Data
    # Check for primary amine/amide
    if not candidate["has_primary_amine"]:
        return f"Incorrect. The IR spectrum shows two N-H bands, characteristic of a primary amine (-NH2). Candidate {llm_answer} ({candidate['name']}) has a secondary amide, which would only show one N-H band."
    
    # Check carbonyl type
    if candidate["carbonyl_type"] == "amide" and ir_data["C=O_stretch"] > 1700:
        return f"Incorrect. The IR peak at 1720 cm-1 is characteristic of an ester. Candidate {llm_answer} ({candidate['name']}) is an amide, which typically has a C=O stretch below 1700 cm-1."

    # 3. Check NMR Data
    # Check for ethoxy group (-O-CH2-CH3) vs other ethyl types
    if not candidate["has_ethoxy_group"]:
        return f"Incorrect. The NMR shows a quartet at a high chemical shift (4.5 ppm), indicating an ethyl group attached to an oxygen (-O-CH2-CH3). Candidate {llm_answer} ({candidate['name']}) does not have this feature."
    
    # Check aromatic substitution pattern
    if candidate["aromatic_substitution"] != nmr_data["aromatic_pattern"]:
        return f"Incorrect. The NMR spectrum shows two doublets in the aromatic region, which is a classic pattern for a para-disubstituted ring. Candidate {llm_answer} ({candidate['name']}) has {candidate['aromatic_substitution']} substitution, which would produce a more complex pattern."

    # If all checks for the selected candidate pass, the answer is correct.
    # Let's double-check the logic for candidate A.
    # Formula: C9H11NO2 -> Correct
    # IR: Primary amine (yes), Ester C=O at 1720 (yes) -> Correct
    # NMR: Ethoxy group (yes), para-substitution (yes) -> Correct
    # All data points are consistent with ethyl 4-aminobenzoate.
    
    return "Correct"

# Run the check
result = check_spectroscopy_answer()
print(result)