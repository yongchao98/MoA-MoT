import math

def estimate_fatigue_life():
    """
    Estimates fatigue life based on the Basquin relationship and Palmgren-Miner rule.
    """
    # Given parameters
    b = -0.09
    
    # Stress level fractions (f) and their corresponding life fractions (n)
    f1, n1 = 1.0, 0.70 # for stress level sigma_e
    f2, n2 = 1.1, 0.20 # for stress level 1.1 * sigma_e
    f3, n3 = 1.2, 0.10 # for stress level 1.2 * sigma_e

    # --- Step 1: Calculate N_e, the life at the endurance limit sigma_e ---
    # We need to assume a ratio for sigma_e / sigma_uts to get a numerical answer.
    # A common rule of thumb for steel is sigma_e / sigma_uts = 0.5.
    sigma_e_over_uts = 0.5
    
    # The fatigue life N for a given stress sigma_a is derived from the two points:
    # (sigma_uts, N=0.5) and (sigma_e, N=N_e)
    # The general formula is N = 0.5 * (sigma_a / sigma_uts)^(1/b)
    # So, for N_e (where sigma_a = sigma_e):
    # N_e = 0.5 * (sigma_e / sigma_uts)^(1/b)
    
    exponent = 1.0 / b
    n_e = 0.5 * (sigma_e_over_uts)**exponent
    
    print("--- Calculation of N_total ---")
    print("Formula: N_total = N_e / (n1 * (f1)^(-1/b) + n2 * (f2)^(-1/b) + n3 * (f3)^(-1/b))")
    print("\nStep 1: Calculate N_e (life at endurance limit)")
    print(f"Assuming sigma_e / sigma_uts = {sigma_e_over_uts}")
    print(f"N_e = 0.5 * ({sigma_e_over_uts})^(1 / {b}) = {n_e:.4f}")
    
    # --- Step 2: Calculate the denominator of the Palmgren-Miner damage sum ---
    # This represents the total damage in one "block" of life fractions.
    # Denominator = n1/((f1)^(1/b)) + n2/((f2)^(1/b)) + n3/((f3)^(1/b))
    # This simplifies to: n1*(f1)^(-1/b) + n2*(f2)^(-1/b) + n3*(f3)^(-1/b)

    term1 = n1 * (f1)**(-exponent)
    term2 = n2 * (f2)**(-exponent)
    term3 = n3 * (f3)**(-exponent)
    
    denominator = term1 + term2 + term3
    
    print("\nStep 2: Calculate the denominator terms based on Palmgren-Miner rule")
    print(f"Term 1 (for {n1*100}% life at {f1}*sigma_e): {n1:.2f} * ({f1})^(-1/{b}) = {term1:.4f}")
    print(f"Term 2 (for {n2*100}% life at {f2}*sigma_e): {n2:.2f} * ({f2})^(-1/{b}) = {term2:.4f}")
    print(f"Term 3 (for {n3*100}% life at {f3}*sigma_e): {n3:.2f} * ({f3})^(-1/{b}) = {term3:.4f}")
    print(f"Denominator = {term1:.4f} + {term2:.4f} + {term3:.4f} = {denominator:.4f}")

    # --- Step 3: Calculate the total fatigue life N_total ---
    n_total = n_e / denominator

    print("\nStep 3: Calculate total fatigue life N_total")
    print(f"N_total = N_e / Denominator = {n_e:.4f} / {denominator:.4f}")
    print(f"N_total = {n_total:.4f} cycles")
    
    # --- Step 4: Format the output to Aeb scientific notation with 3 decimal places ---
    exponent_10 = math.floor(math.log10(n_total))
    mantissa = n_total / (10**exponent_10)
    
    # The notation 'Aeb' implies 'A * 10^b'
    formatted_answer = f"{mantissa:.3f}e{exponent_10}"
    print("\n--- Final Answer ---")
    print(f"The estimated fatigue life in the format Aeb is: {formatted_answer}")

    return formatted_answer

# Execute the function and capture the final answer for the '<<<>>>' block
final_answer = estimate_fatigue_life()

# Final output block as required by instructions
print(f"<<<{final_answer}>>>")
