import math

def estimate_fatigue_life():
    """
    Estimates fatigue life using the Palmgren-Miner rule and Basquin's law.
    """
    # --- Step 1: Define the given parameters and assumptions ---
    
    # Basquin exponent
    b = -0.09
    
    # Life fractions for each stress level
    f1 = 0.70  # at sigma_e
    f2 = 0.20  # at 1.1 * sigma_e
    f3 = 0.10  # at 1.2 * sigma_e
    
    # Assumption: The fatigue endurance limit corresponds to a life (N_e) of 10^7 cycles.
    # This is a standard engineering convention for materials without a distinct fatigue knee.
    N_e = 10**7
    
    # --- Step 2: Set up the Palmgren-Miner equation to solve for total life L ---
    # The equation is: L = N_e / (f1 + f2*(1.1)**(-1/b) + f3*(1.2)**(-1/b))
    
    # Calculate the exponent term
    neg_inv_b = -1.0 / b
    
    # --- Step 3: Calculate each term in the denominator ---
    term1 = f1
    term2_multiplier = math.pow(1.1, neg_inv_b)
    term3_multiplier = math.pow(1.2, neg_inv_b)
    
    term2 = f2 * term2_multiplier
    term3 = f3 * term3_multiplier
    
    denominator = term1 + term2 + term3
    
    # --- Step 4: Calculate the total fatigue life L ---
    L = N_e / denominator
    
    # --- Step 5: Print the results in a clear format ---
    print("Fatigue Life Estimation using Palmgren-Miner Rule")
    print("-" * 50)
    print(f"Basquin exponent (b): {b}")
    print(f"Assumed cycles to failure at endurance limit (N_e): {N_e:e}")
    print("\nPalmgren-Miner Damage Equation for Total Life (L):")
    print(f"L = N_e / (f1 + f2 * (1.1)^(-1/b) + f3 * (1.2)^(-1/b))")
    
    # Print the equation with all numbers plugged in
    print("\nSubstituting the numerical values:")
    final_equation = f"L = {N_e:.0f} / ({f1} + {f2} * (1.1)^({neg_inv_b:.3f}) + {f3} * (1.2)^({neg_inv_b:.3f}))"
    print(final_equation)

    print("\nCalculating the denominator terms:")
    print(f"Term 2: {f2} * {term2_multiplier:.3f} = {term2:.3f}")
    print(f"Term 3: {f3} * {term3_multiplier:.3f} = {term3:.3f}")
    print(f"Total Denominator = {term1} + {term2:.3f} + {term3:.3f} = {denominator:.3f}")

    print("\n--- Final Calculation ---")
    print(f"Total Estimated Fatigue Life (L) = {N_e:.0f} / {denominator:.3f}")
    print(f"L = {L:.3f} cycles")
    
    # Format the final answer as requested: Aeb with 3 decimal places
    final_answer_formatted = f"{L:.3e}"
    print(f"\nFinal Answer (in 'Aeb' format): {final_answer_formatted}")
    
    return final_answer_formatted

# Execute the function and capture the final answer for the required format
final_answer = estimate_fatigue_life()

# The final answer in the required format for submission
# <<<f"{final_answer}">>>