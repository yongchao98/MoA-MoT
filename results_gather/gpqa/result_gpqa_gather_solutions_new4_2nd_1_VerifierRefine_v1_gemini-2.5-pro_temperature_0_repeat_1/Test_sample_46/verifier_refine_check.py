import re

class CompoundIdentifier:
    def __init__(self, question_data, options, llm_answer):
        self.data = question_data
        self.options = options
        self.llm_answer_text = llm_answer
        self.analysis_log = []

    def get_llm_choice(self):
        match = re.search(r'<<<([A-D])>>>', self.llm_answer_text)
        if match:
            return match.group(1)
        return None

    def check_dou(self):
        # C9H11NO2
        C, H, N = 9, 11, 1
        dou = C + 1 - (H / 2) + (N / 2)
        expected_dou = 5
        if dou == expected_dou:
            self.analysis_log.append(f"PASS: Degree of Unsaturation is {dou}, consistent with a benzene ring (4) and a carbonyl (1).")
            return True
        else:
            self.analysis_log.append(f"FAIL: Calculated DoU is {dou}, expected {expected_dou}.")
            return False

    def check_ir(self, candidate):
        features = candidate['features']
        ir_peaks = self.data['ir']
        
        # Check for primary amine
        nh2_peaks = [p for p in ir_peaks if 3300 <= p <= 3500]
        if 'primary_amine' in features:
            if len(nh2_peaks) == 2:
                self.analysis_log.append(f"  - IR (N-H): PASS. Found 2 peaks in 3300-3500 cm-1 range, consistent with a primary amine.")
            else:
                self.analysis_log.append(f"  - IR (N-H): FAIL. Expected 2 N-H peaks for primary amine, but data has {len(nh2_peaks)}.")
                return False
        elif 'secondary_amide' in features:
            if len(nh2_peaks) == 1:
                 self.analysis_log.append(f"  - IR (N-H): PASS. Found 1 peak in 3300-3500 cm-1 range, consistent with a secondary amide.")
            else:
                self.analysis_log.append(f"  - IR (N-H): FAIL. Expected 1 N-H peak for secondary amide, but data has {len(nh2_peaks)}.")
                return False
        else: # e.g. primary amide
             if len(nh2_peaks) != 2:
                self.analysis_log.append(f"  - IR (N-H): FAIL. Expected 2 N-H peaks for primary amide, but data has {len(nh2_peaks)}.")
                return False

        # Check for carbonyl group
        co_peak = next((p for p in ir_peaks if 1650 <= p <= 1800), None)
        if not co_peak:
            self.analysis_log.append(f"  - IR (C=O): FAIL. No carbonyl peak found in the 1650-1800 cm-1 range.")
            return False

        if 'conjugated_ester' in features:
            if 1715 <= co_peak <= 1730:
                self.analysis_log.append(f"  - IR (C=O): PASS. Peak at {co_peak} cm-1 is consistent with a conjugated ester.")
            else:
                self.analysis_log.append(f"  - IR (C=O): FAIL. Peak at {co_peak} cm-1 is outside the expected range (1715-1730) for a conjugated ester.")
                return False
        elif 'amide' in features:
            if 1650 <= co_peak <= 1690:
                self.analysis_log.append(f"  - IR (C=O): PASS. Peak at {co_peak} cm-1 is consistent with an amide.")
            else:
                self.analysis_log.append(f"  - IR (C=O): FAIL. Peak at {co_peak} cm-1 is outside the expected range (1650-1690) for an amide. Data suggests an ester.")
                return False
        elif 'aryl_ester' in features: # e.g. 4-aminophenyl propionate
            if 1750 <= co_peak <= 1770:
                 self.analysis_log.append(f"  - IR (C=O): PASS. Peak at {co_peak} cm-1 is consistent with an aryl ester.")
            else:
                self.analysis_log.append(f"  - IR (C=O): FAIL. Peak at {co_peak} cm-1 is outside the expected range (1750-1770) for an aryl ester. Data suggests a conjugated ester.")
                return False
        
        return True

    def check_nmr(self, candidate):
        features = candidate['features']
        nmr = self.data['nmr']

        # Check substitution pattern
        aromatic_doublets = [s for s, (mult, h) in nmr.items() if 6.5 <= s <= 8.5 and mult == 'd' and h == 2]
        if features['substitution'] == 'para':
            if len(aromatic_doublets) == 2:
                self.analysis_log.append(f"  - NMR (Aromatic): PASS. Two doublets found, consistent with para-substitution.")
            else:
                self.analysis_log.append(f"  - NMR (Aromatic): FAIL. Expected two doublets for para-substitution.")
                return False
        elif features['substitution'] == 'meta':
            if len(aromatic_doublets) != 2:
                self.analysis_log.append(f"  - NMR (Aromatic): PASS. Aromatic pattern is not two doublets, consistent with meta-substitution.")
            else:
                self.analysis_log.append(f"  - NMR (Aromatic): FAIL. Expected complex pattern for meta-substitution, but found two doublets.")
                return False

        # Check for primary amine protons
        amine_proton = any(s == 4.0 and mult == 'bs' and h == 2 for s, (mult, h) in nmr.items())
        if 'primary_amine' in features:
            if amine_proton:
                self.analysis_log.append(f"  - NMR (Amine): PASS. Broad singlet at 4.0 ppm (2H) is consistent with -NH2.")
            else:
                self.analysis_log.append(f"  - NMR (Amine): FAIL. Expected signal for -NH2 not found.")
                return False
        
        # Check alkyl chain
        quartet = next((s for s, (mult, h) in nmr.items() if mult == 'q' and h == 2), None)
        if not quartet:
            self.analysis_log.append(f"  - NMR (Alkyl): FAIL. No quartet found.")
            return False

        if features['alkyl_chain'] == 'ethyl_ester': # -COOCH2CH3
            if 4.1 <= quartet <= 4.6:
                self.analysis_log.append(f"  - NMR (Alkyl): PASS. Quartet at {quartet} ppm is consistent with an ethyl ester (-O-CH2-).")
            else:
                self.analysis_log.append(f"  - NMR (Alkyl): FAIL. Quartet at {quartet} ppm is not in the expected range (4.1-4.6) for an ethyl ester.")
                return False
        elif features['alkyl_chain'] == 'propionate_ester': # -OCOCH2CH3
            if 2.2 <= quartet <= 2.6:
                self.analysis_log.append(f"  - NMR (Alkyl): PASS. Quartet at {quartet} ppm is consistent with a propionate ester (-CO-CH2-).")
            else:
                self.analysis_log.append(f"  - NMR (Alkyl): FAIL. Quartet at {quartet} ppm is not in the expected range (2.2-2.6) for a propionate ester. Data shows {quartet} ppm.")
                return False
        elif features['alkyl_chain'] == 'ethoxy_group': # -OCH2CH3
            if 3.9 <= quartet <= 4.2:
                self.analysis_log.append(f"  - NMR (Alkyl): PASS. Quartet at {quartet} ppm is consistent with an ethoxy group.")
            else:
                self.analysis_log.append(f"  - NMR (Alkyl): FAIL. Quartet at {quartet} ppm is not in the expected range (3.9-4.2) for an ethoxy group. Data shows {quartet} ppm.")
                return False

        return True

    def run_analysis(self):
        llm_choice = self.get_llm_choice()
        if not llm_choice:
            return "Could not parse the answer choice from the LLM's response."

        self.analysis_log.append("Starting analysis...")
        self.check_dou()
        
        correct_candidate_key = None
        for key, candidate in self.options.items():
            self.analysis_log.append(f"\n--- Checking Option {key}: {candidate['name']} ---")
            # All candidates have the correct formula C9H11NO2
            ir_ok = self.check_ir(candidate)
            nmr_ok = self.check_nmr(candidate)
            
            if ir_ok and nmr_ok:
                self.analysis_log.append(f"CONCLUSION: Option {key} is fully consistent with all data.")
                if correct_candidate_key:
                    # This case should not happen with this data, but it's good practice
                    return "Error: Multiple candidates seem to fit the data. Re-check constraints."
                correct_candidate_key = key
            else:
                self.analysis_log.append(f"CONCLUSION: Option {key} is inconsistent with the data.")

        if not correct_candidate_key:
            return "Error: No candidate perfectly matches the provided spectral data."

        if llm_choice == correct_candidate_key:
            return "Correct"
        else:
            llm_candidate_name = self.options[llm_choice]['name']
            correct_candidate_name = self.options[correct_candidate_key]['name']
            
            # Find the specific failure reason for the LLM's choice
            failure_reasons = []
            for log_entry in self.analysis_log:
                if f"Checking Option {llm_choice}" in log_entry:
                    start_index = self.analysis_log.index(log_entry)
                    for entry in self.analysis_log[start_index:]:
                        if "FAIL" in entry:
                            failure_reasons.append(entry.strip())
                        if f"CONCLUSION: Option {llm_choice}" in entry:
                            break
            
            reason_str = ". ".join(failure_reasons)
            return (f"Incorrect. The provided answer was {llm_choice} ({llm_candidate_name}), but the correct answer is {correct_candidate_key} ({correct_candidate_name}).\n"
                    f"Reasoning: The LLM's choice {llm_choice} is incorrect because it fails the following checks: {reason_str}")


# --- Main execution ---
def check_answer():
    # Data from the question
    question_data = {
        'formula': 'C9H11NO2',
        'ir': [3420, 3325, 1720],
        'nmr': {1.20: ('t', 3), 4.0: ('bs', 2), 4.5: ('q', 2), 7.0: ('d', 2), 8.0: ('d', 2)}
    }

    # Options from the question with their expected features
    options = {
        'A': {'name': 'N-(4-ethoxyphenyl)formamide', 'features': {'formula': 'C9H11NO2', 'secondary_amide': True, 'amide': True, 'substitution': 'para', 'alkyl_chain': 'ethoxy_group'}},
        'B': {'name': '3-ethoxybenzamide', 'features': {'formula': 'C9H11NO2', 'primary_amine': True, 'amide': True, 'substitution': 'meta', 'alkyl_chain': 'ethoxy_group'}},
        'C': {'name': '4-aminophenyl propionate', 'features': {'formula': 'C9H11NO2', 'primary_amine': True, 'aryl_ester': True, 'substitution': 'para', 'alkyl_chain': 'propionate_ester'}},
        'D': {'name': 'ethyl 4-aminobenzoate', 'features': {'formula': 'C9H11NO2', 'primary_amine': True, 'conjugated_ester': True, 'substitution': 'para', 'alkyl_chain': 'ethyl_ester'}}
    }

    # The final answer provided by the LLM
    llm_answer = """
    Here is a step-by-step analysis to identify the compound using the provided spectroscopic data.

    **1. Analysis of the Molecular Formula and Degree of Unsaturation (DoU)**
    *   **Molecular Formula**: C9H11NO2
    *   **DoU Calculation**: The degree of unsaturation is calculated using the formula: `DoU = C + 1 - (H/2) + (N/2)`.
        *   DoU = 9 + 1 - (11/2) + (1/2) = 10 - 5.5 + 0.5 = 5.
    *   **Careful Point**: A DoU of 5 is a strong indication of a benzene ring (which accounts for 4 degrees of unsaturation: one ring and three pi bonds) plus one additional pi bond. Given the two oxygen atoms, this is most likely a carbonyl group (C=O).

    **2. Analysis of the IR Spectrum**
    *   **3420 cm⁻¹ and 3325 cm⁻¹**: The presence of two distinct, medium-to-strong bands in this region is the classic signature of a **primary amine (-NH₂)**. These two bands correspond to the asymmetric and symmetric N-H stretching vibrations. This finding makes options with secondary amines or amides unlikely.
    *   **1720 cm⁻¹**: This is a strong absorption in the carbonyl (C=O) region.
    *   **Careful Point**: The frequency of 1720 cm⁻¹ is characteristic of an **ester** that is conjugated with an aromatic ring. It is generally too high for a conjugated amide (which typically appears below 1700 cm⁻¹).
    *   **Conclusion from IR**: The compound contains a primary amine (-NH₂) and a conjugated ester group.

    **3. Analysis of the ¹H NMR Spectrum**
    *   **7.0 ppm (d, 2H) and 8.0 ppm (d, 2H)**:
        *   These signals are in the aromatic region. The pattern of two doublets, each integrating to two protons, is the unmistakable signature of a **1,4-disubstituted (para) benzene ring**.
        *   **Careful Point**: The different chemical shifts indicate that one substituent is electron-donating (shielding the adjacent protons, causing them to appear upfield at 7.0 ppm) and the other is electron-withdrawing (deshielding the adjacent protons, causing them to appear downfield at 8.0 ppm).
    *   **1.20 ppm (t, 3H) and 4.5 ppm (q, 2H)**:
        *   This triplet-quartet combination is the classic pattern for an **ethyl group (-CH₂CH₃)**.
        *   **Careful Point**: The chemical shift of the quartet at **4.5 ppm** is very important. This highly downfield shift for a -CH₂- group indicates it is directly attached to a very electronegative atom, specifically the oxygen of an ester group **(-O-CH₂CH₃)**.
    *   **4.0 ppm (bs, 2H)**:
        *   A broad singlet integrating to two protons is characteristic of exchangeable protons. Combined with the IR data, this signal confirms the presence of the **primary amine (-NH₂) group**.

    **4. Assembling the Structure and Evaluating the Options**
    *   **Summary of Fragments**:
        1.  A 1,4-disubstituted benzene ring.
        2.  A primary amine group (-NH₂), which is electron-donating.
        3.  An ethyl ester group (-COOCH₂CH₃), which is electron-withdrawing.
    *   **Assembled Structure**: Placing these fragments together leads to **ethyl 4-aminobenzoate**. The electron-donating -NH₂ group and the electron-withdrawing -COOCH₂CH₃ group are para to each other, which perfectly explains the NMR aromatic pattern and the chemical shifts.
    *   **Evaluating the Options from the Prompt**:
        *   A) N-(4-ethoxyphenyl)formamide: Incorrect. This is a secondary amide, which contradicts the IR data for a primary amine.
        *   B) 3-ethoxybenzamide: Incorrect. This is an amide (wrong IR C=O frequency) and is 1,3-disubstituted (wrong NMR aromatic pattern).
        *   C) 4-aminophenyl propionate: Incorrect. In this isomer (Ar-O-CO-CH₂CH₃), the NMR signal for the -CH₂- group would be around 2.5 ppm, not 4.5 ppm.
        *   D) ethyl 4-aminobenzoate: **Correct**. This structure perfectly matches all the data: the formula, the primary amine, the ethyl ester, the 1,4-disubstitution, and the electronic effects on the aromatic protons.

    <<<D>>>
    """

    checker = CompoundIdentifier(question_data, options, llm_answer)
    result = checker.run_analysis()
    print(result)

check_answer()