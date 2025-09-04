import re

def check_answer_correctness(llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given spectroscopy problem.

    The function analyzes the spectral data provided in the question and evaluates
    if the chosen candidate molecule matches all the constraints.

    Args:
        llm_answer_text (str): The full text of the LLM's response, which includes
                               the final answer in the format <<<X>>>.

    Returns:
        str: "Correct" if the answer is correct, or a string explaining the
             reason for the incorrectness.
    """

    # --- Step 1: Define constraints from the spectral data ---
    # IR: 3420 & 3325 cm-1 -> Two N-H stretches indicate a primary amine (-NH2).
    # IR: 1720 cm-1 -> Strong C=O stretch, characteristic of a conjugated ester.
    # 1H NMR: 1.20(t,3H) & 4.5(q,2H) -> Ethyl group attached to an ester oxygen (-O-CH2CH3).
    #         The quartet at 4.5 ppm is key, as it's too downfield for -CO-CH2-.
    # 1H NMR: 4.0(bs,2H) -> Broad singlet for 2 protons confirms the primary amine (-NH2).
    # 1H NMR: 7.0(d,2H) & 8.0(d,2H) -> Classic pattern for a para (1,4) disubstituted benzene ring.

    # --- Step 2: Define properties of the candidate molecules ---
    candidates = {
        'A': {
            'name': '3-ethoxybenzamide',
            'is_correct': False,
            'reason': "Incorrect. Candidate A (3-ethoxybenzamide) is meta-substituted, but the NMR shows a para-substitution pattern (two doublets). Also, the IR C=O stretch at 1720 cm-1 is too high for an amide (typically ~1680 cm-1)."
        },
        'B': {
            'name': '4-aminophenyl propionate',
            'is_correct': False,
            'reason': "Incorrect. Candidate B (4-aminophenyl propionate) has a propionyl group (-OCOCH2CH3). Its NMR would show a quartet for the -CH2- group around 2.5 ppm (next to the carbonyl), not at 4.5 ppm as observed."
        },
        'C': {
            'name': 'N-(4-ethoxyphenyl)formamide',
            'is_correct': False,
            'reason': "Incorrect. Candidate C (N-(4-ethoxyphenyl)formamide) is a secondary amide (-NH-), but the IR spectrum clearly shows two N-H stretches, which is characteristic of a primary amine (-NH2)."
        },
        'D': {
            'name': 'ethyl 4-aminobenzoate',
            'is_correct': True,
            'reason': "This candidate is consistent with all data: primary amine (IR), ethyl ester (IR, NMR), and para-substitution (NMR)."
        }
    }
    
    # --- Step 3: Extract the LLM's choice and check it ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find a valid answer in the format <<<X>>> in the provided text."

    chosen_option = match.group(1)

    if candidates[chosen_option]['is_correct']:
        return "Correct"
    else:
        return candidates[chosen_option]['reason']

# The provided answer from the LLM
llm_answer = """
Here is a step-by-step analysis to determine the correct structure for the compound C9H11NO2 based on the provided spectroscopic data.

**1. Analysis of the Molecular Formula and Degree of Unsaturation (DoU)**
*   **Formula:** C9H11NO2
*   **DoU Calculation:** The degree of unsaturation is calculated as `C + 1 - (H/2) + (N/2)`.
    *   DoU = 9 + 1 - (11/2) + (1/2) = 10 - 5.5 + 0.5 = 5.
*   **Careful Point:** A DoU of 5 is a strong indicator of a benzene ring (which accounts for 4 degrees of unsaturation) plus one additional pi bond (e.g., a C=O carbonyl group).

**2. Analysis of the IR Spectrum**
*   **3420 cm⁻¹ and 3325 cm⁻¹:** The presence of two distinct, medium-to-strong bands in this region is the classic signature of a primary amine (-NH₂), corresponding to its asymmetric and symmetric N-H stretching vibrations.
*   **1720 cm⁻¹:** This strong band is in the carbonyl (C=O) stretching region. The frequency is characteristic of an ester that is conjugated with an aromatic ring. It is generally too high for an amide (~1680 cm⁻¹) and fits perfectly for this type of ester.
*   **Careful Point:** The IR data strongly suggests the presence of both a primary amine (-NH₂) and a conjugated ester functional group.

**3. Analysis of the ¹H NMR Spectrum**
*   **1.20 ppm (t, 3H) and 4.5 ppm (q, 2H):** This combination of a triplet (3H) and a quartet (2H) is the definitive pattern for an ethyl group (-CH₂CH₃).
    *   **Careful Point:** The chemical shift of the quartet at 4.5 ppm is very downfield. This indicates that the methylene (-CH₂-) group is attached to a highly electronegative atom, specifically the oxygen of an ester (-O-CH₂CH₃).
*   **4.0 ppm (bs, 2H):** A broad singlet integrating to two protons is characteristic of the exchangeable protons of a primary amine (-NH₂). This confirms the finding from the IR spectrum.
*   **7.0 ppm (d, 2H) and 8.0 ppm (d, 2H):** Two doublets in the aromatic region, each integrating to two protons, is the classic pattern for a 1,4-disubstituted (para) benzene ring.
    *   **Careful Point:** The difference in chemical shifts indicates one substituent is electron-donating (shielding the adjacent protons, causing the upfield signal at 7.0 ppm) and the other is electron-withdrawing (deshielding the adjacent protons, causing the downfield signal at 8.0 ppm). The -NH₂ group is a strong electron-donator, and the ester group (-COOEt) is an electron-withdrawer. This arrangement perfectly explains the observed signals.

**4. Assembling the Structure and Evaluating the Options**
The combined data points to a structure with:
*   A 1,4-disubstituted (para) benzene ring.
*   A primary amine (-NH₂) group.
*   An ethyl ester (-COOCH₂CH₃) group.

This assembly leads to the structure of **ethyl 4-aminobenzoate**.

Let's check the given options:
*   A) 3-ethoxybenzamide: Incorrect. This is a 1,3- (meta) substituted amide. The NMR aromatic pattern and IR C=O stretch would be different.
*   B) 4-aminophenyl propionate: Incorrect. This is an isomer where the ethyl group is part of a propionate attached to a phenolic oxygen. The NMR for the -CH₂- next to the C=O would be a quartet around 2.5 ppm, not 4.5 ppm.
*   C) N-(4-ethoxyphenyl)formamide: Incorrect. This is a secondary amide, which contradicts the IR data (would show only one N-H stretch) and NMR data (would show only one NH proton and a formyl proton).
*   D) ethyl 4-aminobenzoate: Correct. This structure is perfectly consistent with the molecular formula, IR spectrum, and ¹H NMR spectrum.

<<<D>>>
"""

# Run the check
result = check_answer_correctness(llm_answer)
print(result)