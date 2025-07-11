def analyze_drift_barrier(population_name, Ne, s):
    """
    This function demonstrates the drift-barrier hypothesis.

    It calculates the strength of genetic drift (approximated as 1/Ne) and compares it
    to the strength of selection (s) for a beneficial allele (e.g., one that
    lowers the mutation rate).

    Args:
        population_name (str): The name of the example population.
        Ne (float): The effective population size.
        s (float): The selection coefficient (the fitness advantage of the allele).
    """
    strength_of_drift = 1 / Ne
    
    print(f"--- Analysis for {population_name} ---")
    print(f"Effective Population Size (Ne): {Ne:,.0f}")
    print(f"Selective Advantage (s): {s}")
    print(f"Strength of Genetic Drift (1/Ne): {strength_of_drift:.10f}")
    
    if s > strength_of_drift:
        print("Result: Selection is stronger than drift. The beneficial allele is likely to be fixed.")
    else:
        print("Result: Genetic drift is stronger than or equal to selection. The allele's fate is dominated by drift, not selection.")
    print("-" * 50)

# Scenario 1: A large population, like bacteria.
# Here, even very weak selection is effective.
Ne_bacteria = 1e8  # 100,000,000
s_weak = 1e-7      # 0.0000001 (a very small advantage)
analyze_drift_barrier("Bacterial Population", Ne_bacteria, s_weak)

# Scenario 2: A small population, like a vertebrate.
# Here, drift is much stronger and can overwhelm the same weak selective force.
Ne_vertebrate = 1e4 # 10,000
analyze_drift_barrier("Vertebrate Population", Ne_vertebrate, s_weak)
