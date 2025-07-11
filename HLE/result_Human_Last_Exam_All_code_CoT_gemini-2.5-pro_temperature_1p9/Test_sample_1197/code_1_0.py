import math

def illustrate_selection_challenge():
    """
    This script illustrates why a non-linear relationship between genome architecture
    and selection pressure challenges simple evolutionary models.

    Concept:
    The fate of a slightly deleterious mutation (e.g., an increase in intron length)
    depends on the balance between selection (s) and genetic drift. The effectiveness
    of purifying selection is often estimated by the product Ne*s, where Ne is the
    effective population size.
    - If Ne*s << 1, drift dominates, and the mutation behaves as if neutral.
    - If Ne*s >> 1, selection dominates, and the mutation is purged.

    The challenge (from Answer D):
    If the selection coefficient 's' is not constant, but instead a non-linear function
    of genomic context (like gene density), then simple models assuming a constant 's'
    will fail to predict evolutionary outcomes across the genome.

    This code compares the effectiveness of selection (Ne*s) for a linear vs. a
    non-linear model of 's' as a function of gene density.
    """
    
    # --- Parameters ---
    # Effective population size
    Ne = 5000
    # A base selection coefficient (small and negative, as it's deleterious)
    s_base = -0.00001
    
    print("--- Illustrating the Challenge to Drift Models ---")
    print(f"Assuming Effective Population Size (Ne) = {Ne}")
    print(f"Base selection coefficient against long introns (s_base) = {s_base}\n")
    print("Comparing Linear vs. Non-Linear effect of Gene Density on Selection\n")
    
    print(f"{'Gene Density':<15} | {'Linear s':<15} | {'Linear Ne*s':<15} | {'Non-Linear s':<15} | {'Non-Linear Ne*s':<15}")
    print("-" * 80)
    
    # Iterate through different levels of gene density (from 0 to 1.0)
    for i in range(11):
        density = i / 10.0
        
        # --- Model 1: Linear Relationship ---
        # Selection 's' increases linearly with gene density.
        s_linear = s_base * density
        selection_effectiveness_linear = Ne * abs(s_linear) # abs() because we care about magnitude
        
        # --- Model 2: Non-Linear Relationship ---
        # Selection 's' increases as the fourth power of density. A sharp, non-linear increase.
        s_nonlinear = s_base * (density ** 4)
        selection_effectiveness_nonlinear = Ne * abs(s_nonlinear)

        print(f"{density:<15.2f} | {s_linear:<15.6f} | {selection_effectiveness_linear:<15.2f} | {s_nonlinear:<15.6f} | {selection_effectiveness_nonlinear:<15.2f}")
        
    print("-" * 80)
    print("\nConclusion:")
    print("In the non-linear model, the effectiveness of selection skyrockets at high gene densities.")
    print("A simple model assuming a linear relationship would fail to predict this,")
    print("thus demonstrating the challenge posed by such complex, non-linear interactions.")


if __name__ == "__main__":
    illustrate_selection_challenge()
