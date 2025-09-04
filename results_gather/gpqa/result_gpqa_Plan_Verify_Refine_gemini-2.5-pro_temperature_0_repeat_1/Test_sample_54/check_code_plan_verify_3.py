import re
from collections import Counter

def check_correctness_of_nmr_answer():
    """
    This function checks the correctness of the LLM's answer for an NMR spectroscopy problem.
    It codifies the rules of NMR interpretation for the given options and compares them
    against the provided spectral data.
    """

    # --- Part 1: Define the problem statement and the LLM's answer ---

    # The 1H NMR data from the question
    given_data = {
        'signals': [
            {'ppm': 7.0, 'integration': 1, 'splitting': 'd', 'J': 16.0},
            {'ppm': 5.5, 'integration': 1, 'splitting': 'dq'},
            {'ppm': 2.1, 'integration': 3, 'splitting': 's'},
            {'ppm': 1.6, 'integration': 3, 'splitting': 'd'}
        ],
        'options': {
            'A': 'Trans-propenyl acetate',
            'B': 'Cis-propenyl acetate',
            'C': 'Cis-butenyl acetate',
            'D': 'Trans-butenyl acetate'
        }
    }

    # The LLM's full response, from which we'll extract the final answer
    llm_response = """
**Plan:**
1.  Analyze the coupling constant (J value) of the vinyl proton signal at 7.0 ppm to determine the stereochemistry (cis or trans) of the double bond.
2.  Analyze the integration and splitting patterns of the remaining signals to distinguish between a propenyl and a butenyl group.
3.  Combine the findings to identify the correct compound from the given options.

**Action:**
The most diagnostic signal for determining stereochemistry is the one at 7.0 ppm, which is a doublet with a coupling constant of J = 16.0 Hz. This J value arises from the coupling between the two protons on the double bond. I will compare this value to the standard ranges for cis and trans alkenes.
<search>J coupling constant ranges for cis and trans alkenes</search>
<information>
- **Trans-coupling (³J_trans):** The coupling between two protons in a trans configuration on a double bond is typically in the range of **11-18 Hz**.
- **Cis-coupling (³J_cis):** The coupling between two protons in a cis configuration on a double bond is typically in the range of **6-15 Hz**, most commonly **6-12 Hz**.
</information>
**Plan:**
1.  The J-coupling constant of 16.0 Hz confirms a *trans* configuration, which eliminates the *cis* options.
2.  Analyze the remaining signals (5.5 ppm, 2.1 ppm, 1.6 ppm) to differentiate between the trans-propenyl and trans-butenyl structures.
3.  Conclude the final structure.

**Action:**
The observed coupling constant of J = 16.0 Hz falls squarely in the typical range for **trans** protons on a double bond (11-18 Hz). This allows us to confidently eliminate the *cis* isomers, options B and C.

Now we must distinguish between the two remaining options: A) Trans-propenyl acetate and D) Trans-butenyl acetate. Let's analyze the rest of the spectrum:

*   **2.1 ppm (3H, s):** A 3-proton singlet is characteristic of the methyl group of an acetate (`CH₃-C=O`), which is present in both remaining options.
*   **1.6 ppm (3H, d):** A 3-proton doublet indicates a methyl group (`CH₃-`) coupled to a single adjacent proton. This pattern is consistent with a propenyl group (`CH₃-CH=`). In contrast, a butenyl group (`CH₃-CH₂-CH=`) would have a terminal methyl group appearing as a triplet (coupled to the adjacent `-CH₂-` group). This signal confirms the presence of a propenyl group.
*   **5.5 ppm (1H, dq):** A 1-proton doublet of quartets indicates a proton coupled to one other proton (giving the doublet splitting) and to a group of three protons (giving the quartet splitting). This perfectly describes the vinyl proton in the `CH₃-CH=` fragment. It is split into a doublet by the other vinyl proton and into a quartet by the three protons of the methyl group.

The combination of a *trans* coupling constant (16.0 Hz), a 3H doublet (1.6 ppm), and a 1H doublet of quartets (5.5 ppm) is uniquely consistent with the structure of Trans-propenyl acetate.

<<<A>>>
    """

    # --- Part 2: Define the "knowledge base" for the candidate compounds ---

    knowledge_base = {
        'Trans-propenyl acetate': {
            'total_protons': 8,
            'j_coupling_range': (11, 18), # Trans
            'expected_patterns': [
                {'integration': 3, 'splitting': 's'},  # Acetate CH3
                {'integration': 3, 'splitting': 'd'},  # Vinyl CH3
                {'integration': 1, 'splitting': 'dq'}, # Vinyl H next to CH3
                {'integration': 1, 'splitting': 'd'}   # Vinyl H next to O
            ]
        },
        'Cis-propenyl acetate': {
            'total_protons': 8,
            'j_coupling_range': (6, 15), # Cis
            'expected_patterns': [
                {'integration': 3, 'splitting': 's'},
                {'integration': 3, 'splitting': 'd'},
                {'integration': 1, 'splitting': 'dq'},
                {'integration': 1, 'splitting': 'd'}
            ]
        },
        'Trans-butenyl acetate': {
            'total_protons': 10, # Different proton count
            'j_coupling_range': (11, 18), # Trans
            'expected_patterns': [ # Different patterns
                {'integration': 3, 'splitting': 's'},  # Acetate CH3
                {'integration': 3, 'splitting': 't'},  # Terminal CH3
                {'integration': 2, 'splitting': 'm'},  # Methylene CH2
                {'integration': 1, 'splitting': 'm'},  # Vinyl H 1
                {'integration': 1, 'splitting': 'd'}   # Vinyl H 2
            ]
        },
        'Cis-butenyl acetate': {
            'total_protons': 10,
            'j_coupling_range': (6, 15), # Cis
            'expected_patterns': [
                {'integration': 3, 'splitting': 's'},
                {'integration': 3, 'splitting': 't'},
                {'integration': 2, 'splitting': 'm'},
                {'integration': 1, 'splitting': 'm'},
                {'integration': 1, 'splitting': 'd'}
            ]
        }
    }

    # --- Part 3: Extract and validate the LLM's chosen answer ---

    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Failure: Could not find a final answer in the format <<<X>>>."
    
    answer_key = match.group(1)
    chosen_compound_name = given_data['options'].get(answer_key)
    if not chosen_compound_name:
        return f"Failure: The extracted answer key '{answer_key}' is not a valid option."

    # --- Part 4: Check the chosen answer against the data and knowledge base ---

    expected = knowledge_base[chosen_compound_name]

    # Check 1: Total proton count
    given_total_protons = sum(s['integration'] for s in given_data['signals'])
    if given_total_protons != expected['total_protons']:
        return (f"Incorrect. The answer is {chosen_compound_name}, which should have {expected['total_protons']} protons. "
                f"However, the sum of integrations in the NMR data is {given_total_protons}H.")

    # Check 2: J-coupling constant for stereochemistry
    j_value_signal = next((s for s in given_data['signals'] if 'J' in s), None)
    j_value = j_value_signal['J']
    j_min, j_max = expected['j_coupling_range']
    if not (j_min <= j_value <= j_max):
        stereochemistry = 'trans' if 'Trans' in chosen_compound_name else 'cis'
        return (f"Incorrect. The answer is {chosen_compound_name}. This implies a {stereochemistry} configuration, "
                f"which typically has a J-coupling constant in the range of {j_min}-{j_max} Hz. "
                f"The given J-value of {j_value} Hz is outside this range.")

    # Check 3: Signal integration and splitting patterns
    given_patterns = [{'integration': s['integration'], 'splitting': s['splitting']} for s in given_data['signals']]
    expected_counts = Counter(tuple(sorted(d.items())) for d in expected['expected_patterns'])
    given_counts = Counter(tuple(sorted(d.items())) for d in given_patterns)

    if expected_counts != given_counts:
        return (f"Incorrect. The integration and splitting patterns in the data do not match those expected for {chosen_compound_name}. "
                f"Expected patterns: {expected['expected_patterns']}. "
                f"Given patterns: {given_patterns}.")

    # --- Part 5: Final Verdict ---
    return "Correct"

# Execute the check and print the result
result = check_correctness_of_nmr_answer()
print(result)