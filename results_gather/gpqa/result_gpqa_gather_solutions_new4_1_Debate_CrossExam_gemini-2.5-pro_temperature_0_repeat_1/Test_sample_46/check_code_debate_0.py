import re

def check_molecular_formula(name, expected_formula="C9H11NO2"):
    """Calculates the molecular formula from a chemical name and checks it."""
    formulas = {
        "ethyl 4-aminobenzoate": "C9H11NO2",
        "4-aminophenyl propionate": "C9H11NO2",
        "3-ethoxybenzamide": "C9H11NO2",
        "N-(4-ethoxyphenyl)formamide": "C9H11NO2"
    }
    actual_formula = formulas.get(name)
    if actual_formula != expected_formula:
        return f"Molecular formula mismatch for {name}. Expected {expected_formula}, but it is {actual_formula}."
    return None

def check_ir_spectrum(features, ir_data):
    """Checks if the compound's features match the IR data."""
    reasons = []
    
    # Check for primary amine vs secondary amide
    has_primary_amine_peaks = 3300 < ir_data['nh_peak1'] < 3500 and 3300 < ir_data['nh_peak2'] < 3500
    if features['amine_type'] == 'primary':
        if not has_primary_amine_peaks:
            reasons.append("Expected two N-H stretch peaks for a primary amine in the 3300-3500 cm-1 range, which were not specified as present.")
    elif features['amine_type'] == 'secondary_amide':
        if has_primary_amine_peaks:
            reasons.append("Data shows two N-H stretch peaks (3420, 3325 cm-1), characteristic of a primary amine, not a secondary amide which has only one.")
            
    # Check carbonyl peak
    co_peak = ir_data['co_peak']
    expected_co_range = features['co_range']
    if not (expected_co_range[0] <= co_peak <= expected_co_range[1]):
        reasons.append(f"Observed C=O peak at {co_peak} cm-1 is outside the expected range of {expected_co_range[0]}-{expected_co_range[1]} cm-1 for a {features['carbonyl_type']}.")
        
    return reasons

def check_nmr_spectrum(features, nmr_data):
    """Checks if the compound's features match the 1H NMR data."""
    reasons = []
    
    # Check aromatic substitution pattern
    is_para_pattern = ('d' in nmr_data['aromatic1_mult'] and 'd' in nmr_data['aromatic2_mult'] and
                       nmr_data['aromatic1_h'] == 2 and nmr_data['aromatic2_h'] == 2)
    if features['substitution'] == '1,4-para':
        if not is_para_pattern:
            reasons.append("Expected a 1,4-para substitution pattern (two doublets, 2H each) in the aromatic region, which was not observed.")
    elif features['substitution'] == '1,3-meta':
        if is_para_pattern:
            reasons.append("Observed a 1,4-para substitution pattern (two doublets), which is inconsistent with the expected complex pattern for a 1,3-meta substituted ring.")

    # Check ethyl group type
    quartet_shift = nmr_data['quartet_ppm']
    if features['ethyl_type'] == 'ethoxy': # -O-CH2-CH3
        if not (4.0 <= quartet_shift <= 4.6):
            reasons.append(f"Observed quartet at {quartet_shift} ppm is not in the expected range (4.0-4.6 ppm) for an ethoxy group (-O-CH2-).")
    elif features['ethyl_type'] == 'propionyl': # -C(=O)-CH2-CH3
        if not (2.2 <= quartet_shift <= 2.8):
             reasons.append(f"Observed quartet at {quartet_shift} ppm is not in the expected range (2.2-2.8 ppm) for a propionyl group (-C(=O)-CH2-). The observed shift points to an ethoxy group.")

    # Check amine/amide protons
    amine_protons = nmr_data['broad_singlet_h']
    if features['amine_type'] == 'primary' and amine_protons != 2:
        reasons.append(f"Expected a signal for 2 amine protons (-NH2), but observed a signal for {amine_protons}H.")
    elif features['amine_type'] == 'secondary_amide' and amine_protons != 1:
        # Note: The question data has a 2H broad singlet, so this check is relevant.
        reasons.append(f"Expected a signal for 1 amide proton (-NH-), but the data shows a broad singlet for 2H.")
        
    return reasons

def check_correctness():
    """
    Main function to check the correctness of the answer by analyzing all candidate structures.
    """
    # --- Data from the question ---
    question_data = {
        "formula": "C9H11NO2",
        "ir": {
            "nh_peak1": 3420,
            "nh_peak2": 3325,
            "co_peak": 1720
        },
        "nmr": {
            "triplet_ppm": 1.20, "triplet_h": 3,
            "broad_singlet_ppm": 4.0, "broad_singlet_h": 2,
            "quartet_ppm": 4.5, "quartet_h": 2,
            "aromatic1_ppm": 7.0, "aromatic1_mult": "d", "aromatic1_h": 2,
            "aromatic2_ppm": 8.0, "aromatic2_mult": "d", "aromatic2_h": 2,
        }
    }
    
    # --- Expected features of candidate structures ---
    candidates = {
        "A": {
            "name": "ethyl 4-aminobenzoate",
            "features": {
                "amine_type": "primary",
                "carbonyl_type": "conjugated ester",
                "co_range": (1715, 1730),
                "substitution": "1,4-para",
                "ethyl_type": "ethoxy" # -O-CH2-CH3
            }
        },
        "B": {
            "name": "4-aminophenyl propionate",
            "features": {
                "amine_type": "primary",
                "carbonyl_type": "phenyl ester",
                "co_range": (1730, 1750), # Phenyl esters are slightly higher
                "substitution": "1,4-para",
                "ethyl_type": "propionyl" # -C(=O)-CH2-CH3
            }
        },
        "C": {
            "name": "3-ethoxybenzamide",
            "features": {
                "amine_type": "primary_amide", # -CONH2 has 2 N-H peaks
                "carbonyl_type": "primary amide",
                "co_range": (1650, 1690),
                "substitution": "1,3-meta",
                "ethyl_type": "ethoxy"
            }
        },
        "D": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "features": {
                "amine_type": "secondary_amide", # -NH- has 1 N-H peak
                "carbonyl_type": "secondary amide",
                "co_range": (1650, 1690),
                "substitution": "1,4-para",
                "ethyl_type": "ethoxy"
            }
        }
    }
    
    llm_answer = "A"
    
    correct_candidate = None
    all_reasons = {}

    for key, data in candidates.items():
        name = data["name"]
        features = data["features"]
        
        errors = []
        
        # Check 1: Molecular Formula
        formula_error = check_molecular_formula(name, question_data["formula"])
        if formula_error:
            errors.append(formula_error)
        
        # Check 2: IR Spectrum
        ir_errors = check_ir_spectrum(features, question_data["ir"])
        errors.extend(ir_errors)
        
        # Check 3: NMR Spectrum
        nmr_errors = check_nmr_spectrum(features, question_data["nmr"])
        errors.extend(nmr_errors)
        
        if not errors:
            correct_candidate = key
        all_reasons[key] = errors

    if correct_candidate is None:
        return "The provided data is inconsistent with all options, or there is an error in the analysis logic."

    if llm_answer == correct_candidate:
        return "Correct"
    else:
        reasoning = f"The provided answer '{llm_answer}' is incorrect. The correct answer is '{correct_candidate}'.\n"
        reasoning += f"Here is why '{llm_answer}' ({candidates[llm_answer]['name']}) is wrong:\n"
        for reason in all_reasons[llm_answer]:
            reasoning += f"- {reason}\n"
        reasoning += f"\nHere is why '{correct_candidate}' ({candidates[correct_candidate]['name']}) is correct:\n"
        reasoning += "- It is a primary amine, matching the two IR N-H peaks at 3420 and 3325 cm-1.\n"
        reasoning += "- It is a conjugated ester, matching the strong IR C=O peak at 1720 cm-1.\n"
        reasoning += "- It has a 1,4-para substituted ring, matching the two doublets (2H each) in the NMR aromatic region.\n"
        reasoning += "- It has an ethyl ester group (-O-CH2CH3), matching the NMR quartet at 4.5 ppm and triplet at 1.20 ppm.\n"
        reasoning += "- It has a primary amine (-NH2), matching the broad singlet for 2H at 4.0 ppm."
        return reasoning

# Run the check
print(check_correctness())