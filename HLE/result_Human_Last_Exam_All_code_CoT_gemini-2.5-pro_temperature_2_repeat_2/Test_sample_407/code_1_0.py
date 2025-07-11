import numpy as np

def calculate_gene_flow_cost():
    """
    Simulates the analysis of a yeast gene flow experiment.

    This function uses hypothetical data to demonstrate how to calculate the cost
    of gene flow by measuring the selection coefficient (s) of F1 hybrids
    and an F2 generation (produced by "within mating") relative to a
    parental ("no gene flow") line.
    """
    # --- Step 1: Define hypothetical experimental data ---
    # These values represent a fitness proxy, like growth rate (doublings/hour).
    parental_growth_rate = 0.50
    f1_hybrid_growth_rate = 0.48  # F1 hybrids may show a slight fitness cost
    f2_generation_growth_rate = 0.42 # F2 may show a larger cost due to outbreeding depression

    print("--- Hypothetical Experimental Data ---")
    print(f"Parental Line (No Gene Flow) Growth Rate: {parental_growth_rate}")
    print(f"F1 Hybrid Growth Rate: {f1_hybrid_growth_rate}")
    print(f"F2 Generation (from Hybrid within-mating) Growth Rate: {f2_generation_growth_rate}\n")

    # --- Step 2: Calculate the Selection Coefficient of the F1 Hybrids ---
    # The selection coefficient (s) measures the fitness cost relative to the reference population.
    # Relative fitness (W) = (Fitness of genotype) / (Fitness of reference)
    # Selection coefficient (s) = 1 - W
    
    print("--- Analysis of F1 Hybrids (Immediate Cost) ---")
    # Calculate F1 relative fitness
    W_f1 = f1_hybrid_growth_rate / parental_growth_rate
    # Calculate F1 selection coefficient
    s_f1 = 1 - W_f1
    
    print("1. Calculate relative fitness (W) of the F1 hybrid:")
    print(f"   W_f1 = F1 Hybrid Growth Rate / Parental Growth Rate")
    print(f"   W_f1 = {f1_hybrid_growth_rate} / {parental_growth_rate} = {W_f1:.2f}")

    print("\n2. Calculate the selection coefficient (s) against the F1 hybrid:")
    print(f"   s_f1 = 1 - W_f1")
    print(f"   s_f1 = 1 - {W_f1:.2f} = {s_f1:.2f}")
    print(f"The immediate cost of gene flow (in F1) is a selection coefficient of {s_f1:.2f}.\n")


    # --- Step 3: Account for the effects of meiosis by analyzing the F2 generation ---
    # This reveals the cost due to incompatible gene combinations after recombination.

    print("--- Analysis of F2 Generation (Cost after Meiosis) ---")
    # Calculate F2 relative fitness
    W_f2 = f2_generation_growth_rate / parental_growth_rate
    # Calculate F2 selection coefficient
    s_f2 = 1 - W_f2

    print("1. Calculate relative fitness (W) of the F2 generation:")
    print(f"   W_f2 = F2 Generation Growth Rate / Parental Growth Rate")
    print(f"   W_f2 = {f2_generation_growth_rate} / {parental_growth_rate} = {W_f2:.2f}")

    print("\n2. Calculate the selection coefficient (s) against the F2 generation:")
    print(f"   s_f2 = 1 - W_f2")
    print(f"   s_f2 = 1 - {W_f2:.2f} = {s_f2:.2f}")
    print(f"The cost after meiosis (outbreeding depression) results in a selection coefficient of {s_f2:.2f}.")

calculate_gene_flow_cost()