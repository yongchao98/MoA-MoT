import math

def check_answer():
    """
    Checks the correctness of the answer to the astrophysics redshift question.
    """
    # --- Define constants and given values ---

    # Rest wavelength of Lyman-alpha line in Angstroms
    lambda_rest = 1216  # Å

    # Atmospheric cutoff for ground-based optical telescopes.
    # Light with wavelengths shorter than this is absorbed by the atmosphere (mainly ozone).
    # A typical value used for this limit is between 3200 Å and 3500 Å.
    # We will use the value from the explanation to directly check its logic.
    lambda_cutoff = 3500  # Å

    # The options provided in the question
    options = {
        "A": 3.0,
        "B": 1.9,
        "C": 2.4,
        "D": 1.2
    }
    
    # The answer given by the LLM
    llm_answer_key = "B"

    # --- Verification Logic ---

    # The question asks for the "lower limit", so we should find the smallest
    # redshift 'z' from the options that makes the Lyman-alpha line observable.
    # We sort the options by their redshift value to check them in increasing order.
    sorted_options = sorted(options.items(), key=lambda item: item[1])

    correct_key = None
    
    print("Verifying the lower limit for redshift...")
    print(f"Rest Wavelength (Lyman-alpha): {lambda_rest} Å")
    print(f"Assumed Atmospheric Cutoff Wavelength: {lambda_cutoff} Å\n")

    for key, z in sorted_options:
        # Calculate the observed wavelength using the redshift formula
        lambda_observed = lambda_rest * (1 + z)
        
        print(f"Checking option {key} (z = {z}):")
        print(f"  λ_observed = {lambda_rest} * (1 + {z}) = {lambda_observed:.1f} Å")

        # Check if the observed wavelength is greater than the atmospheric cutoff
        if lambda_observed >= lambda_cutoff:
            print(f"  Result: {lambda_observed:.1f} Å >= {lambda_cutoff} Å. This is DETECTABLE.")
            # Since we are checking in increasing order of z, the first one we find
            # that is detectable is the correct lower limit among the choices.
            if correct_key is None:
                correct_key = key
                print(f"  -> This is the first detectable option, so it's the lower limit.")
        else:
            print(f"  Result: {lambda_observed:.1f} Å < {lambda_cutoff} Å. This is NOT DETECTABLE from the ground.")
        print("-" * 20)

    # --- Final Conclusion ---
    if correct_key is None:
        return "Error in checking logic or problem constraints. No suitable answer found."

    print(f"\nBased on the calculation, the correct option is '{correct_key}'.")
    print(f"The LLM's answer was '{llm_answer_key}'.")

    if llm_answer_key == correct_key:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer was {llm_answer_key}, but the calculation shows "
                f"that {correct_key} is the lowest redshift option for which the Lyman-alpha line "
                f"({lambda_rest} Å) is shifted past the atmospheric cutoff of ~{lambda_cutoff} Å.")

# Run the check
result = check_answer()
print(f"\nFinal Verdict: {result}")