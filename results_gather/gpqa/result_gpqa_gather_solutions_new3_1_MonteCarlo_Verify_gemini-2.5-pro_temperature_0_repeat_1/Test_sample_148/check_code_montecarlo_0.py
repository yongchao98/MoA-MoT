import re

def check_answer_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer based on the problem's constraints.

    The function simulates the logical deduction process:
    1.  Analyzes the mass spectrometry data to determine if the species are isomers.
    2.  Analyzes the NMR and LC data to determine if the species are distinguishable.
    3.  Compares these findings with the properties of enantiomers, diastereoisomers,
        precursors, and side-products to find the only matching explanation.
    """

    # --- Step 1: Define the experimental evidence from the question ---
    evidence = {
        "are_isomers": True,  # From MS: both peaks have the same mass as the expected molecule.
        "are_distinguishable_by_nmr": True,  # From NMR: two distinct peaks for one proton.
        "are_separable_by_lc": True  # From LC: two clearly defined peaks.
    }

    # --- Step 2: Define the properties of each possible explanation (the options) ---
    # Note: 'Standard' or 'achiral' conditions are assumed for NMR and LC.
    options_properties = {
        'A': {  # Enantiomers
            "name": "a mixture of enantiomers",
            "are_isomers": True,
            "are_distinguishable_by_nmr": False,
            "are_separable_by_lc": False
        },
        'B': {  # Diastereoisomers
            "name": "a mixture of diastereoisomers",
            "are_isomers": True,
            "are_distinguishable_by_nmr": True,
            "are_separable_by_lc": True
        },
        'C': {  # 'Double coupling' product
            "name": "a 'double coupling' product",
            "are_isomers": False,  # This would be a different molecule with a higher mass.
        },
        'D': {  # Precursor contamination
            "name": "contamination with a precursor",
            "are_isomers": False,  # This would be a different molecule with a lower mass.
        }
    }

    # --- Step 3: Extract the chosen option from the LLM's answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not find a valid answer in the format <<<X>>> in the provided text."
    
    chosen_option_key = match.group(1)
    
    # --- Step 4: Determine the correct answer based on the evidence ---
    correct_option_key = None
    for key, properties in options_properties.items():
        # Check if the properties of the option match the evidence
        if (properties.get("are_isomers") == evidence["are_isomers"] and
            properties.get("are_distinguishable_by_nmr") == evidence["are_distinguishable_by_nmr"] and
            properties.get("are_separable_by_lc") == evidence["are_separable_by_lc"]):
            correct_option_key = key
            break

    # --- Step 5: Compare the LLM's choice with the correct answer and provide feedback ---
    if chosen_option_key == correct_option_key:
        return "Correct"
    else:
        chosen_properties = options_properties[chosen_option_key]
        
        # Provide a specific reason for the error
        if not chosen_properties.get("are_isomers", False):
            return (f"Incorrect. The chosen answer ({chosen_option_key}) suggests {chosen_properties['name']}, which is not an isomer of the target molecule. "
                    "This contradicts the mass spectrometry data, which shows both species have the same mass as the expected product, proving they are isomers.")
        
        if not chosen_properties.get("are_distinguishable_by_nmr", False):
            return (f"Incorrect. The chosen answer ({chosen_option_key}) suggests {chosen_properties['name']}. Under standard (achiral) conditions, these are not distinguishable by NMR. "
                    "This contradicts the NMR data, which shows two distinct peaks for the alpha-proton.")

        if not chosen_properties.get("are_separable_by_lc", False):
            return (f"Incorrect. The chosen answer ({chosen_option_key}) suggests {chosen_properties['name']}. Under standard (achiral) conditions, these are not separable by LC. "
                    "This contradicts the LC data, which shows two distinct, separable peaks.")

        return f"Incorrect. The chosen option {chosen_option_key} does not align with all experimental evidence."


# The final answer from the LLM response is passed to the checking function.
llm_final_answer = "<<<B>>>"
result = check_answer_correctness(llm_final_answer)
print(result)