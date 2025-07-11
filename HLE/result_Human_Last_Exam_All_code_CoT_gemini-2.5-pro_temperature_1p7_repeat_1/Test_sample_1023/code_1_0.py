def calculate_extinction_rate_factor():
    """
    Calculates how much greater the extinction rate for a morphospecies is
    compared to the extinction rate for an evolutionary species based on the problem's assumptions.
    """

    # For calculation purposes, let's assign a baseline value to the true extinction rate `mu`.
    # The actual value is irrelevant as it will cancel out in the final ratio.
    mu = 1.0

    # --- Step 1: Define rates for Evolutionary Species ---
    # The extinction rate for an evolutionary species is the true extinction rate.
    mu_e = mu
    
    print(f"Let's assume the true extinction rate of an evolutionary species (μ_e) is {mu_e}")

    # --- Step 2: Apply assumptions to find relationships between rates ---

    # Assumption 1: Steady state of names implies the true speciation rate equals the true extinction rate.
    # lambda = mu
    lmbda = mu
    print(f"From the steady-state assumption, the true speciation rate (λ) is equal to μ, so λ = {lmbda}")

    # Assumption 2: Symmetry in pseudoextinction implies the rate of anagenesis equals the rate of bifurcation.
    # The rate of pseudoextinction from bifurcation is 0.5 * lambda.
    # So, lambda_anagenesis = 0.5 * lambda
    lambda_anagenesis = 0.5 * lmbda
    print(f"From the symmetry assumption, the anagenesis rate (λ_a) is 0.5 * λ, so λ_a = {lambda_anagenesis}")

    # --- Step 3: Calculate the extinction rate for Morphospecies ---

    # The extinction rate for a morphospecies (μ_m) is the sum of:
    # 1. True extinction rate (μ)
    # 2. Pseudoextinction from bifurcation (0.5 * λ)
    # 3. Pseudoextinction from anagenesis (λ_a)
    
    mu_component = mu
    bifurcation_component = 0.5 * lmbda
    anagenesis_component = lambda_anagenesis

    mu_m = mu_component + bifurcation_component + anagenesis_component

    print("\nThe extinction rate for a morphospecies (μ_m) is the sum of:")
    print(f" - True extinction (μ): {mu_component}")
    print(f" - Pseudoextinction from bifurcation (0.5 * λ): {bifurcation_component}")
    print(f" - Pseudoextinction from anagenesis (λ_a): {anagenesis_component}")
    
    # --- Step 4: Calculate and print the final factor ---
    
    factor = mu_m / mu_e
    
    print("\nSo, the total extinction rate for a morphospecies is:")
    print(f"μ_m = {mu_component} + {bifurcation_component} + {anagenesis_component} = {mu_m}")
    
    print("\nThe final ratio of the extinction rates is:")
    print(f"Ratio = μ_m / μ_e = {mu_m} / {mu_e} = {factor}")
    
    return factor

# Run the calculation and store the final answer
final_answer = calculate_extinction_rate_factor()
print(f"\n<<<The extinction rate for a morphospecies is {final_answer} times greater than for an evolutionary species.>>>")
print(f"\n<<<{final_answer}>>>")