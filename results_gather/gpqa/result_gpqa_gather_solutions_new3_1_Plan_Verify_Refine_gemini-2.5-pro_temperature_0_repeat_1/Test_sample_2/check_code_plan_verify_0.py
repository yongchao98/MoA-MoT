import numpy as np
import re

def check_quantum_expectation_value():
    """
    This function calculates the expectation value for the given quantum mechanics problem
    and checks if the provided LLM answer is correct.
    """
    # Define the coefficients of the quantum state |ψ⟩ = c_up|↑⟩ + c_down|↓⟩
    c_up = 0.5
    c_down = np.sqrt(3) / 2

    # --- Step 1: Verify state normalization ---
    # The sum of the squared magnitudes of the coefficients should be 1.
    # This is a good sanity check.
    normalization = np.abs(c_up)**2 + np.abs(c_down)**2
    if not np.isclose(normalization, 1.0):
        return f"State is not normalized. |c_up|^2 + |c_down|^2 = {normalization}"

    # --- Step 2: Calculate individual expectation values ---
    # For a state |ψ⟩ = c_up|↑⟩ + c_down|↓⟩:
    # ⟨σ_z⟩ = |c_up|^2 - |c_down|^2
    # ⟨σ_x⟩ = 2 * Re(c_up* * c_down)
    
    exp_sigma_z = np.abs(c_up)**2 - np.abs(c_down)**2
    
    # Since coefficients are real, Re(c_up* * c_down) simplifies to c_up * c_down
    exp_sigma_x = 2 * c_up * c_down

    # --- Step 3: Calculate the total expectation value ---
    # The operator is O = 10*σ_z + 5*σ_x
    # ⟨O⟩ = 10*⟨σ_z⟩ + 5*⟨σ_x⟩
    total_exp_exact = 10 * exp_sigma_z + 5 * exp_sigma_x

    # --- Step 4: Round to one decimal place as per the question ---
    calculated_value_rounded = round(total_exp_exact, 1)

    # --- Step 5: Compare with the LLM's answer ---
    # The options given in the question are:
    # A) 0.85, B) -0.7, C) -1.4, D) 1.65
    options = {'A': 0.85, 'B': -0.7, 'C': -1.4, 'D': 1.65}
    
    # The final answer from the LLM to be checked
    llm_final_answer_text = "<<<B>>>" # Extracted from the provided response

    # Extract the letter from the LLM's answer
    match = re.search(r'<<<([A-D])>>>', llm_final_answer_text)
    if not match:
        return "Could not parse the final answer from the LLM's response."
        
    llm_choice_letter = match.group(1)
    llm_choice_value = options.get(llm_choice_letter)

    # Check if the LLM's choice corresponds to the calculated value
    if np.isclose(calculated_value_rounded, llm_choice_value):
        return "Correct"
    else:
        # Find the correct letter choice based on our calculation
        correct_letter = None
        for letter, value in options.items():
            if np.isclose(calculated_value_rounded, value):
                correct_letter = letter
                break
        
        reason = (f"Incorrect. The calculated expectation value is {total_exp_exact:.4f}, "
                  f"which rounds to {calculated_value_rounded}. "
                  f"This corresponds to option {correct_letter} ({options[correct_letter]}). "
                  f"The LLM answered with option {llm_choice_letter} ({llm_choice_value}).")
        return reason

# Execute the check and print the result
result = check_quantum_expectation_value()
print(result)