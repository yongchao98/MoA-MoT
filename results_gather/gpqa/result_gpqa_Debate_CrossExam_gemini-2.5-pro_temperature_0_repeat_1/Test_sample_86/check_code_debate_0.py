import re

def check_correctness_of_llm_answer():
    """
    This function checks the correctness of the LLM's answer by verifying its reasoning
    against the constraints provided in the chemistry question.
    """

    # --- Step 1: Define constraints from the question ---
    # - Compound is a di-substituted 6-membered aromatic ring.
    # - Total of 8 carbon atoms.
    # - FTIR indicates a carbonyl group (C=O).
    # - FTIR indicates an aromatic-halogen bond (Ar-X).
    # - This implies the structure is X-C6H4-R, where R contains the remaining atoms.
    
    # --- Step 2: Verify the proposed structure from the LLM's reasoning ---
    # The LLM proposes a halo-acetophenone structure: X-C6H4-COCH3.
    # Let's check this structure against the constraints.
    # - Carbons: 6 (ring) + 1 (carbonyl C) + 1 (methyl C) = 8. This matches.
    # - Di-substituted aromatic ring: Yes, a halogen (X) and an acetyl group (-COCH3). This matches.
    # - Carbonyl group: Yes, the acetyl group has one. This matches.
    # - Aromatic-halogen bond: Yes, the halogen is directly on the ring. This matches.
    # The proposed structure is valid. The LLM's reasoning is sound.

    # --- Step 3: Analyze the NMR data for the chosen answer (C) ---
    llm_choice = 'C'
    options_data = {
        'A': "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        'B': "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        'C': "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        'D': "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)"
    }
    
    # A simple parser for the NMR data string
    def parse_nmr(data_str):
        signals = re.findall(r'(\d+\.\d+)\s*\((\d+)H,\s*([a-z]+)\)', data_str)
        return [{'ppm': float(s[0]), 'H': int(s[1]), 'split': s[2]} for s in signals]

    # --- Step 4: Verify if the NMR data for option C matches the proposed structure ---
    # The expected structure is para-halo-acetophenone.
    # Expected NMR features:
    # 1. Aromatic region (4H): Two doublets, each integrating to 2H.
    # 2. Aliphatic region (3H): One singlet for the acetyl methyl group, around 2.1-2.6 ppm.
    
    c_signals = parse_nmr(options_data['C'])
    
    # Check 1: Aromatic protons
    aromatic_signals = [s for s in c_signals if 6.5 <= s['ppm'] <= 8.5]
    aromatic_proton_count = sum(s['H'] for s in aromatic_signals)
    
    if aromatic_proton_count != 4:
        return f"Reason for incorrectness: The chosen answer C is wrong. A di-substituted benzene ring must have 4 aromatic protons, but the data shows {aromatic_proton_count}."

    # Check 2: Para-substitution pattern
    is_para_pattern = (
        len(aromatic_signals) == 2 and
        all(s['split'] == 'd' for s in aromatic_signals) and
        all(s['H'] == 2 for s in aromatic_signals)
    )
    if not is_para_pattern:
        return "Reason for incorrectness: The chosen answer C is wrong. The aromatic signals do not show the classic para-substitution pattern (two doublets of 2H each) expected from the symmetry."

    # Check 3: Acetyl group protons (-COCH3)
    aliphatic_signals = [s for s in c_signals if s['ppm'] < 6.5]
    
    if len(aliphatic_signals) != 1:
        return f"Reason for incorrectness: The chosen answer C is wrong. It should have only one non-aromatic signal for the acetyl group, but it has {len(aliphatic_signals)}."
        
    acetyl_signal = aliphatic_signals[0]
    is_acetyl_methyl = (
        2.0 <= acetyl_signal['ppm'] <= 2.7 and
        acetyl_signal['H'] == 3 and
        acetyl_signal['split'] == 's'
    )
    
    if not is_acetyl_methyl:
        return f"Reason for incorrectness: The chosen answer C is wrong. The signal at {acetyl_signal['ppm']} ppm does not match an acetyl methyl group (expected: 3H singlet around 2.1-2.6 ppm)."

    # --- Step 5: Briefly verify the rejection of other options to confirm robust reasoning ---
    # Option A: No signals in the aromatic region (6.5-8.5 ppm). Fails.
    a_aromatic_protons = sum(s['H'] for s in parse_nmr(options_data['A']) if 6.5 <= s['ppm'] <= 8.5)
    if a_aromatic_protons != 0:
        return "Code check failed: Logic error in analyzing option A."

    # Option B: Has an aldehyde proton (9.9 ppm). To have 8 carbons and an Ar-X bond, the structure would be p-halo-methylbenzaldehyde. This would have a methyl singlet (~2.5 ppm), not a methylene singlet (3.7 ppm). The LLM's reasoning is sound. Fails.
    b_signals = parse_nmr(options_data['B'])
    is_aldehyde = any(s['ppm'] > 9.0 for s in b_signals)
    has_wrong_aliphatic = any(s['ppm'] == 3.7 and s['H'] == 2 for s in b_signals)
    if not (is_aldehyde and has_wrong_aliphatic):
        return "Code check failed: Logic error in analyzing option B."

    # Option D: Only 1 aromatic proton. A di-substituted ring must have 4. Fails.
    d_aromatic_protons = sum(s['H'] for s in parse_nmr(options_data['D']) if 6.5 <= s['ppm'] <= 8.5)
    if d_aromatic_protons == 4:
        return "Code check failed: Logic error in analyzing option D."

    # --- Final Conclusion ---
    # All checks passed. The LLM correctly identified the structure and matched it to option C,
    # while correctly eliminating the other options.
    return "Correct"

# Execute the check and print the result
result = check_correctness_of_llm_answer()
print(result)