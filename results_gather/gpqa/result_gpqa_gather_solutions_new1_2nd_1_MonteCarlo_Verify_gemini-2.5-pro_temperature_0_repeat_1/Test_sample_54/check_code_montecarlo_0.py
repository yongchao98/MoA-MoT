import re

def check_nmr_answer():
    """
    Checks the correctness of the LLM's answer for the NMR problem.
    """
    # --- Problem Definition ---
    # 1H NMR data provided in the question
    nmr_data = {
        "signals": [
            {"ppm": 7.0, "integration": 1, "multiplicity": "d", "J_Hz": 16.0},
            {"ppm": 5.5, "integration": 1, "multiplicity": "dq"},
            {"ppm": 2.1, "integration": 3, "multiplicity": "s"},
            {"ppm": 1.6, "integration": 3, "multiplicity": "d"},
        ]
    }
    # Options provided in the question
    options = {
        "A": "Trans-butenyl acetate",
        "B": "Cis-butenyl acetate",
        "C": "Cis-propenyl acetate",
        "D": "Trans-propenyl acetate",
    }
    
    # --- Expected Properties of Compounds ---
    # Define the expected properties for each possible compound
    compound_properties = {
        "Cis-propenyl acetate": {"protons": 8, "stereochemistry": "cis", "j_range_hz": (6, 15)},
        "Trans-propenyl acetate": {"protons": 8, "stereochemistry": "trans", "j_range_hz": (11, 18)},
        "Cis-butenyl acetate": {"protons": 10, "stereochemistry": "cis", "j_range_hz": (6, 15)},
        "Trans-butenyl acetate": {"protons": 10, "stereochemistry": "trans", "j_range_hz": (11, 18)},
    }

    # --- LLM's Answer ---
    llm_answer_text = """
    Here is a step-by-step analysis to determine the correct compound from the provided 1H NMR data.

    ### Step-by-step Analysis

    1.  **Analyze the Total Number of Protons (Integration)**
        *   The spectrum shows signals with integrations of 1H, 1H, 3H, and 3H.
        *   The total number of protons in the molecule is 1 + 1 + 3 + 3 = 8 protons.
        *   **Careful Point:** This is a powerful initial filter. We must compare this count to the options.
            *   Propenyl acetate (C₅H₈O₂) has 8 protons. This is a match.
            *   Butenyl acetate (C₆H₁₀O₂) has 10 protons. This is not a match.
        *   **Conclusion:** This analysis immediately eliminates both butenyl acetate options. The unknown compound must be either Cis-propenyl acetate or Trans-propenyl acetate.

    2.  **Analyze the Coupling Constant (J-value) to Determine Stereochemistry**
        *   The signal at 7.0 ppm is a doublet with a coupling constant of **J = 16.0 Hz**. This J-value describes the interaction between the two protons on the carbon-carbon double bond.
        *   **Careful Point:** The magnitude of the coupling constant between vinylic protons is the most definitive piece of information for determining the geometry of the double bond.
            *   A **trans** configuration (protons on opposite sides of the double bond) typically shows a large coupling constant, in the range of **11-18 Hz**.
            *   A **cis** configuration (protons on the same side) shows a smaller coupling constant, typically in the range of **6-15 Hz** (and often below 12 Hz).
        *   **Conclusion:** The observed J-value of **16.0 Hz** falls squarely and unambiguously within the range for a **trans** configuration. This eliminates the cis isomer.

    3.  **Verify the Structure with All Signals**
        *   From the steps above, the only possible structure is **Trans-propenyl acetate**. We can confirm that all signals are consistent with this structure: CH₃-CH=CH-O-C(=O)-CH₃.
        *   **2.1 ppm (3H, s):** A singlet for 3 protons is characteristic of the methyl group of an acetate (CH₃-C=O), which has no adjacent protons. This matches.
        *   **1.6 ppm (3H, d):** A doublet for 3 protons indicates a methyl group (CH₃-) coupled to a single adjacent proton (-CH=). This matches the terminal methyl group of the propenyl chain.
        *   **5.5 ppm (1H, dq):** A doublet of quartets for 1 proton. This proton is split by two different groups: the other vinylic proton (giving the doublet, 'd') and the three protons of the adjacent methyl group (giving the quartet, 'q'). This complex splitting pattern perfectly describes the vinylic proton in the middle of the propenyl chain (=CH-CH₃).
        *   **7.0 ppm (1H, d, J = 16.0 Hz):** The downfield doublet is the vinylic proton attached to the oxygen (-O-CH=), split by the other vinylic proton in a *trans* arrangement. This matches.

    ### Final Conclusion

    The analysis of the NMR data provides two decisive constraints:
    1.  The proton count (8H) limits the options to **propenyl acetate**.
    2.  The large coupling constant (J = 16.0 Hz) confirms the **trans** stereochemistry.

    Both constraints point unequivocally to **Trans-propenyl acetate**. Matching this to the provided options:
    A) Trans-butenyl acetate
    B) Cis-butenyl acetate
    C) Cis-propenyl acetate
    D) Trans-propenyl acetate

    The correct choice is D.

    <<<D>>>
    """

    # 1. Parse the final answer from the text
    match = re.search(r'<<<([A-D])>>>\s*$', llm_answer_text)
    if not match:
        return "Failure: Could not find the final answer in the required format '<<<X>>>' at the end of the response."
    
    llm_choice = match.group(1)
    chosen_compound_name = options.get(llm_choice)

    # --- Verification Logic ---

    # Check 1: Total Proton Count
    data_total_protons = sum(s["integration"] for s in nmr_data["signals"])
    expected_protons = compound_properties[chosen_compound_name]["protons"]
    
    if data_total_protons != expected_protons:
        return (f"Incorrect. The total proton count constraint is not satisfied. "
                f"The NMR data shows {data_total_protons} protons, but the chosen answer "
                f"'{chosen_compound_name}' has {expected_protons} protons.")

    # Check 2: Stereochemistry via J-coupling constant
    vinylic_signal_with_j = next((s for s in nmr_data["signals"] if "J_Hz" in s), None)
    data_j_value = vinylic_signal_with_j["J_Hz"]
    expected_j_range = compound_properties[chosen_compound_name]["j_range_hz"]
    expected_stereochem = compound_properties[chosen_compound_name]["stereochemistry"]

    if not (expected_j_range[0] <= data_j_value <= expected_j_range[1]):
        return (f"Incorrect. The stereochemistry constraint is not satisfied. "
                f"The data's J-coupling constant of {data_j_value} Hz is not in the expected range "
                f"for a '{expected_stereochem}' configuration, which is typically {expected_j_range} Hz.")

    # Check 3: Structural Fragments
    # All answers are acetates, so the 3H singlet is expected.
    # The key difference is propenyl vs. butenyl.
    is_propenyl = "propenyl" in chosen_compound_name.lower()
    if is_propenyl:
        # A propenyl group (CH3-CH=CH-) should have a 3H doublet and a 1H dq.
        has_ch3_doublet = any(s["integration"] == 3 and s["multiplicity"] == 'd' for s in nmr_data["signals"])
        has_ch_dq = any(s["integration"] == 1 and s["multiplicity"] == 'dq' for s in nmr_data["signals"])
        if not (has_ch3_doublet and has_ch_dq):
            return "Incorrect. The splitting patterns in the data do not match a propenyl group."
    
    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_nmr_answer()
print(result)