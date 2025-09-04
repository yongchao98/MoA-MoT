import collections

def check_correctness():
    """
    Checks the correctness of the identified compound based on spectral data.
    The function analyzes each piece of data (Formula, IR, NMR) and verifies
    if the proposed answer (C: ethyl 4-aminobenzoate) uniquely satisfies all constraints.
    """
    # --- Define Spectral Data from the Question ---
    spectral_data = {
        "formula": "C9H11NO2",
        "ir_nh_bands": 2,  # Two bands at 3420 & 3325 cm-1 indicate a primary amine
        "ir_co_freq": 1720,  # Strong band at 1720 cm-1
        "nmr_aromatic_pattern": "para",  # Two doublets (2H each) -> 1,4-disubstitution
        "nmr_amine_protons_shift": 4.0, # bs, 2H -> primary amine protons
        "nmr_ethyl_quartet_shift": 4.5  # q, 2H -> key for ethyl group environment
    }

    # --- Define Properties of Candidate Compounds ---
    # Note: Expected NMR shifts and IR frequencies are approximate ranges.
    candidates = {
        'A': {
            "name": "N-(4-ethoxyphenyl)formamide",
            "formula": "C9H11NO2",
            "amine_type": "secondary_amide", # Expects 1 N-H band
            "carbonyl_type": "amide",
            "carbonyl_freq_range": (1660, 1700),
            "substitution_pattern": "para",
            "ethyl_quartet_shift_ppm": 4.0, # Ethoxy group (-O-CH2-Ar)
            "has_primary_amine_protons": False
        },
        'B': {
            "name": "3-ethoxybenzamide",
            "formula": "C9H11NO2",
            "amine_type": "primary_amide", # Expects 2 N-H bands, but it's an amide
            "carbonyl_type": "amide",
            "carbonyl_freq_range": (1650, 1690),
            "substitution_pattern": "meta",
            "ethyl_quartet_shift_ppm": 4.0, # Ethoxy group (-O-CH2-Ar)
            "has_primary_amine_protons": False
        },
        'C': {
            "name": "ethyl 4-aminobenzoate",
            "formula": "C9H11NO2",
            "amine_type": "primary", # Expects 2 N-H bands
            "carbonyl_type": "ester",
            "carbonyl_freq_range": (1715, 1730), # Conjugated ester
            "substitution_pattern": "para",
            "ethyl_quartet_shift_ppm": 4.3, # Ethyl ester group (-O-CH2-CO)
            "has_primary_amine_protons": True
        },
        'D': {
            "name": "4-aminophenyl propionate",
            "formula": "C9H11NO2",
            "amine_type": "primary", # Expects 2 N-H bands
            "carbonyl_type": "ester",
            "carbonyl_freq_range": (1750, 1770), # Phenyl ester (higher freq)
            "substitution_pattern": "para",
            "ethyl_quartet_shift_ppm": 2.5, # Propionate group (-CO-CH2-CH3)
            "has_primary_amine_protons": True
        }
    }

    # The final answer provided by the LLM to be checked
    llm_answer_id = 'C'
    answer_candidate = candidates[llm_answer_id]
    
    error_messages = []

    # --- Verification Steps ---

    # 1. Check IR N-H stretch (Primary Amine)
    if answer_candidate["amine_type"] != "primary":
        error_messages.append(f"The answer {llm_answer_id} ({answer_candidate['name']}) is a {answer_candidate['amine_type']}, but the IR data (2 N-H bands) indicates a primary amine.")

    # 2. Check IR C=O stretch (Conjugated Ester)
    co_freq = spectral_data['ir_co_freq']
    co_range = answer_candidate['carbonyl_freq_range']
    if not (co_range[0] <= co_freq <= co_range[1]):
        error_messages.append(f"The IR C=O stretch at {co_freq} cm-1 is inconsistent with the expected range for {answer_candidate['name']} ({co_range[0]}-{co_range[1]} cm-1). It is characteristic of a conjugated ester.")

    # 3. Check NMR Aromatic Pattern (Para substitution)
    if answer_candidate["substitution_pattern"] != spectral_data["nmr_aromatic_pattern"]:
        error_messages.append(f"The answer {llm_answer_id} ({answer_candidate['name']}) has {answer_candidate['substitution_pattern']} substitution, but the NMR aromatic pattern (two doublets) indicates para substitution.")

    # 4. Check NMR Ethyl Group Quartet Shift (The most decisive constraint)
    q_shift_observed = spectral_data['nmr_ethyl_quartet_shift']
    q_shift_expected = answer_candidate['ethyl_quartet_shift_ppm']
    # Use a tolerance of +/- 0.3 ppm for this check
    if not (q_shift_expected - 0.3 <= q_shift_observed <= q_shift_expected + 0.3):
        error_messages.append(f"The observed NMR quartet shift of {q_shift_observed} ppm is inconsistent with the structure of {answer_candidate['name']}, which would show a signal around {q_shift_expected} ppm. The observed shift is characteristic of an ethyl group attached to an ester oxygen (-O-CH2-).")

    # 5. Check NMR Primary Amine Protons
    if not answer_candidate["has_primary_amine_protons"]:
        error_messages.append(f"The answer {llm_answer_id} ({answer_candidate['name']}) does not have primary amine protons, but the NMR shows a signal for 2H at {spectral_data['nmr_amine_protons_shift']} ppm.")

    # --- Final Verdict ---
    if error_messages:
        return "Incorrect. The proposed answer 'C' does not satisfy all constraints:\n" + "\n".join(error_messages)
    else:
        # Final check to ensure no other candidate also fits the data
        for other_id, other_candidate in candidates.items():
            if other_id == llm_answer_id:
                continue
            
            passes_all = True
            if other_candidate["amine_type"] != "primary": passes_all = False
            if not (other_candidate['carbonyl_freq_range'][0] <= co_freq <= other_candidate['carbonyl_freq_range'][1]): passes_all = False
            if other_candidate["substitution_pattern"] != spectral_data["nmr_aromatic_pattern"]: passes_all = False
            if not (other_candidate['ethyl_quartet_shift_ppm'] - 0.3 <= q_shift_observed <= other_candidate['ethyl_quartet_shift_ppm'] + 0.3): passes_all = False
            if not other_candidate["has_primary_amine_protons"]: passes_all = False
            
            if passes_all:
                return f"Incorrect. The analysis is ambiguous as both {answer_candidate['name']} ({llm_answer_id}) and {other_candidate['name']} ({other_id}) could satisfy the spectral constraints."

        return "Correct"

# Execute the check and print the result
result = check_correctness()
print(result)