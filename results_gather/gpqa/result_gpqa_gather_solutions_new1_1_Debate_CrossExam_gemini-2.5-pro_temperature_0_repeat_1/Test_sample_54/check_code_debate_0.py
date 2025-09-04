import collections

def check_nmr_answer():
    """
    This function programmatically checks the correctness of the LLM's answer
    by analyzing the provided 1H NMR data against the properties of the candidate compounds.
    """
    # 1. Define the problem data from the question
    nmr_data = {
        "signals": [
            {"ppm": 7.0, "integration": 1, "multiplicity": "d", "J_Hz": 16.0},
            {"ppm": 5.5, "integration": 1, "multiplicity": "dq"},
            {"ppm": 2.1, "integration": 3, "multiplicity": "s"},
            {"ppm": 1.6, "integration": 3, "multiplicity": "d"},
        ],
        "options": {
            "A": "Cis-propenyl acetate",
            "B": "Trans-propenyl acetate",
            "C": "Cis-butenyl acetate",
            "D": "Trans-butenyl acetate",
        }
    }
    
    # The final answer provided by the LLM
    llm_answer = "B"

    # 2. Define the known properties of the candidate compounds
    candidate_properties = {
        "Cis-propenyl acetate": {"protons": 8, "stereochemistry": "cis", "group": "propenyl"},
        "Trans-propenyl acetate": {"protons": 8, "stereochemistry": "trans", "group": "propenyl"},
        "Cis-butenyl acetate": {"protons": 10, "stereochemistry": "cis", "group": "butenyl"},
        "Trans-butenyl acetate": {"protons": 10, "stereochemistry": "trans", "group": "butenyl"},
    }

    # 3. Perform a series of checks based on the NMR data
    
    # Check 1: Total Proton Count
    total_protons_observed = sum(s['integration'] for s in nmr_data['signals'])
    if total_protons_observed != 8:
        return f"Data parsing error: The sum of protons from the NMR data is {total_protons_observed}, not 8."

    possible_candidates = list(nmr_data['options'].values())
    
    candidates_after_proton_check = [
        c for c in possible_candidates if candidate_properties[c]["protons"] == total_protons_observed
    ]
    
    if not candidates_after_proton_check:
        return f"Incorrect: No candidate compound has the observed {total_protons_observed} protons. The butenyl acetates (10H) should be eliminated."

    # Check 2: Stereochemistry from J-coupling constant
    j_coupling_signal = next((s for s in nmr_data['signals'] if s.get("J_Hz")), None)
    if not j_coupling_signal:
        return "Error: No J-coupling constant found in the data to determine stereochemistry."

    j_value = j_coupling_signal['J_Hz']
    
    # Rule: trans coupling is ~12-18 Hz; cis coupling is ~6-12 Hz.
    if 12 <= j_value <= 18:
        observed_stereochemistry = "trans"
    elif 6 <= j_value < 12:
        observed_stereochemistry = "cis"
    else:
        return f"Incorrect: The J-coupling constant of {j_value} Hz is ambiguous or outside typical ranges."

    candidates_after_stereo_check = [
        c for c in candidates_after_proton_check if candidate_properties[c]["stereochemistry"] == observed_stereochemistry
    ]

    if not candidates_after_stereo_check:
        return f"Incorrect: After filtering for {total_protons_observed} protons, no remaining candidate matches the observed '{observed_stereochemistry}' stereochemistry from the J-value of {j_value} Hz."

    # Check 3: Fragmentation pattern (propenyl vs. butenyl)
    # A propenyl group (CH3-CH=) is confirmed by a methyl doublet and a vinyl doublet of quartets.
    has_methyl_doublet = any(s['integration'] == 3 and s['multiplicity'] == 'd' for s in nmr_data['signals'])
    has_vinyl_dq = any(s['integration'] == 1 and s['multiplicity'] == 'dq' for s in nmr_data['signals'])
    
    if has_methyl_doublet and has_vinyl_dq:
        observed_group = "propenyl"
    else:
        # This check is redundant if the proton count is correct, but good for verification
        observed_group = "unknown"

    final_candidates = [
        c for c in candidates_after_stereo_check if candidate_properties[c]["group"] == observed_group
    ]

    # 4. Final conclusion
    if len(final_candidates) != 1:
        return f"Incorrect: The analysis resulted in {len(final_candidates)} possible candidates: {final_candidates}. The data is either ambiguous or inconsistent."

    deduced_compound_name = final_candidates[0]
    
    # Find the letter corresponding to the deduced compound
    deduced_answer_letter = None
    for letter, name in nmr_data['options'].items():
        if name == deduced_compound_name:
            deduced_answer_letter = letter
            break
            
    # 5. Compare with the LLM's answer
    if llm_answer == deduced_answer_letter:
        return "Correct"
    else:
        return (f"Incorrect: The provided answer is {llm_answer} ({nmr_data['options'][llm_answer]}), "
                f"but the correct answer based on the NMR data is {deduced_answer_letter} ({deduced_compound_name}).\n"
                f"Reasoning:\n"
                f"1. The total proton count is {total_protons_observed}H, which eliminates the butenyl acetate options (10H).\n"
                f"2. The J-coupling constant is {j_value} Hz, which is characteristic of a 'trans' double bond, eliminating the 'cis' isomer.\n"
                f"Therefore, the only matching compound is Trans-propenyl acetate.")

# Execute the check and print the result
result = check_nmr_answer()
print(result)