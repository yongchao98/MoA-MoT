import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by systematically evaluating each option
    against the constraints given in the chemistry problem.
    """
    llm_answer = "C"
    options = {
        "A": "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        "B": "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        "C": "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        "D": "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)"
    }

    # Store analysis results for each option
    analysis_results = {}

    for option_letter, nmr_data in options.items():
        # --- Step 1: Parse NMR data string ---
        # Regex to find ppm, integration (H), and multiplicity (s, d, t, etc.)
        signals = re.findall(r"(\d+\.?\d+)\s*\((\d+)H,\s*([a-z]+)\)", nmr_data)
        parsed_signals = [{'ppm': float(s[0]), 'H': int(s[1]), 'mult': s[2]} for s in signals]

        # --- Step 2: Apply Constraints ---
        # Constraint: Di-substituted 6-membered aromatic ring (must have 4 aromatic protons)
        aromatic_protons = sum(s['H'] for s in parsed_signals if 6.5 <= s['ppm'] <= 8.5)
        if aromatic_protons != 4:
            analysis_results[option_letter] = f"Incorrect. A di-substituted benzene ring must have 4 aromatic protons, but this option has {aromatic_protons}."
            continue

        # Constraint: 8 total carbons (C6 from ring + C2 from substituents)
        # Constraint: Carbonyl group (C=O) present
        # Let's identify the non-aromatic part and check its structure
        non_aromatic_signals = [s for s in parsed_signals if not (6.5 <= s['ppm'] <= 8.5)]
        
        is_valid_substituent = False
        proposed_structure = "Unknown"
        shift_plausibility = "Good"

        # Check for acetyl group (-COCH3) -> C2, has C=O
        # Characteristic signal: ~2.1-2.6 ppm, 3H, singlet
        acetyl_signals = [s for s in non_aromatic_signals if 2.1 <= s['ppm'] <= 2.6 and s['H'] == 3 and s['mult'] == 's']
        if len(acetyl_signals) == 1 and len(non_aromatic_signals) == 1:
            is_valid_substituent = True
            proposed_structure = "para-halo-acetophenone"
            # Acetyl is a strong electron-withdrawing group (EWG), justifying downfield aromatic shifts.
            if not all(s['ppm'] > 7.4 for s in parsed_signals if 6.5 <= s['ppm'] <= 8.5):
                 shift_plausibility = "Poor, aromatic shifts should be downfield with an acetyl group."

        # Check for formyl-methyl group (-CH2CHO) -> C2, has C=O
        # Characteristic signals: aldehyde H (~9-10 ppm, 1H, s) and methylene H (~3.5-3.8 ppm, 2H)
        aldehyde_signals = [s for s in non_aromatic_signals if 9.0 <= s['ppm'] <= 10.0 and s['H'] == 1]
        methylene_signals = [s for s in non_aromatic_signals if 3.5 <= s['ppm'] <= 3.8 and s['H'] == 2]
        if len(aldehyde_signals) == 1 and len(methylene_signals) == 1 and len(non_aromatic_signals) == 2:
            is_valid_substituent = True
            proposed_structure = "para-halo-phenylacetaldehyde"
            # -CH2CHO is a weak EWG. Very downfield aromatic shifts (>7.5 ppm) are less likely.
            if all(s['ppm'] > 7.5 for s in parsed_signals if 6.5 <= s['ppm'] <= 8.5):
                 shift_plausibility = "Poor, aromatic shifts are too far downfield for a weak EWG like -CH2CHO."

        if not is_valid_substituent:
            analysis_results[option_letter] = f"Incorrect. The non-aromatic signals do not correspond to a valid C2 substituent containing a carbonyl group."
            continue
        
        if shift_plausibility != "Good":
             analysis_results[option_letter] = f"Plausible but unlikely. Structure is {proposed_structure}. Reason: {shift_plausibility}"
        else:
             analysis_results[option_letter] = f"Correctly fits all criteria. Proposed structure: {proposed_structure}."


    # --- Step 3: Final Evaluation ---
    # Find the best option based on the analysis
    best_option = None
    for letter, result in analysis_results.items():
        if "Correctly fits all criteria" in result:
            best_option = letter
            break # Found the most plausible option

    # If no single best option, check for the most plausible among the remaining
    if not best_option:
        # In this specific problem, C is a perfect fit and B is a poor fit.
        # The logic above correctly identifies C as the best option.
        pass

    if best_option == llm_answer:
        return "Correct"
    else:
        # Provide a detailed reason for the discrepancy
        reason = f"The LLM's answer is {llm_answer}, but the analysis points to {best_option}.\n"
        reason += "Here is the breakdown:\n"
        for letter, result in analysis_results.items():
            reason += f"  - Option {letter}: {result}\n"
        
        # Find the specific reason why the LLM's choice might be wrong or why our analysis differs
        llm_choice_analysis = analysis_results.get(llm_answer, "Analysis not available.")
        if "Incorrect" in llm_choice_analysis or "unlikely" in llm_choice_analysis:
             reason += f"\nThe LLM's choice '{llm_answer}' is likely incorrect because: {llm_choice_analysis}"
        
        return reason.strip()

# Run the check
result = check_correctness()
print(result)