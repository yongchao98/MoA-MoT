import re

def check_nmr_answer():
    """
    Checks the correctness of the LLM's answer by verifying the constraints from the question.
    """
    # --- Step 1: Define constraints from the question ---
    # Formula C8H7XO -> 7 total protons
    # Di-substituted benzene ring -> 4 aromatic protons
    expected_total_protons = 7
    expected_aromatic_protons = 4
    aromatic_region_threshold = 6.5  # ppm

    # --- Step 2: Define the options and the given answer ---
    options = {
        "A": "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        "B": "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        "C": "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        "D": "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)"
    }
    llm_answer = "C"

    # --- Step 3: Analyze the chosen answer 'C' ---
    chosen_spectrum = options.get(llm_answer)
    if not chosen_spectrum:
        return f"Incorrect. The answer '{llm_answer}' is not one of the valid options."

    # Use regex to parse the NMR data string into (shift, integration, multiplicity)
    try:
        peaks = re.findall(r'(\d+\.?\d*)\s*\((\d+)H,\s*([a-z]+)\)', chosen_spectrum)
        peaks = [(float(p[0]), int(p[1]), p[2]) for p in peaks]
    except (ValueError, IndexError):
        return f"Incorrect. Could not parse the NMR data for the chosen answer '{llm_answer}'."

    # Check total proton count
    total_protons = sum(p[1] for p in peaks)
    if total_protons != expected_total_protons:
        return f"Incorrect. The chosen answer '{llm_answer}' has an incorrect total proton count. Expected {expected_total_protons}, but found {total_protons}."

    # Check aromatic proton count
    aromatic_protons = sum(p[1] for p in peaks if p[0] > aromatic_region_threshold)
    if aromatic_protons != expected_aromatic_protons:
        return f"Incorrect. The chosen answer '{llm_answer}' has an incorrect aromatic proton count. Expected {expected_aromatic_protons}, but found {aromatic_protons}."

    # --- Step 4: Check for specific structural features for halo-acetophenone ---
    # The structure should be a halo-acetophenone, which has a -COCH3 group.
    # This should give a 3H singlet around 2.1-2.6 ppm.
    has_methyl_ketone_signal = any(
        p[1] == 3 and p[2] == 's' and 2.0 < p[0] < 2.7 for p in peaks
    )
    if not has_methyl_ketone_signal:
        return f"Incorrect. While the proton counts for answer '{llm_answer}' are correct, the signals do not match the most likely structure (halo-acetophenone). A 3H singlet for the acetyl group is expected but not found in the correct region."

    # --- Step 5: Verify that other options are incorrect ---
    # Option A: 1 aromatic proton, 7 total protons. Fails.
    # Option B: 0 aromatic protons, 7 total protons. Fails.
    # Option D: 4 aromatic protons, 7 total protons. This corresponds to a halo-phenylacetaldehyde.
    # The choice of C over D is a reasonable chemical judgment (acetophenones are very common).
    # Since the chosen answer C fits all constraints and represents the most plausible structure, it is correct.

    return "Correct"

# Run the check
result = check_nmr_answer()
print(result)