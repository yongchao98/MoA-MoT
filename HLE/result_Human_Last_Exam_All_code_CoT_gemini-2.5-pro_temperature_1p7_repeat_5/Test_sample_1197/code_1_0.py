import numpy as np

def genome_model_challenge():
    """
    Simulates why the correlation between synonymous and nonsynonymous substitution rates
    challenges predictive models of genetic drift.
    """
    # --- Step 1: Define Parameters for a "Simple Predictive Model" ---
    # This model assumes a universal mutation rate and predicts the nonsynonymous
    # rate (Ka) based on a constant purifying selection factor.
    # It assumes Ka is independent of the synonymous rate (Ks).
    baseline_mutation_rate = 0.01
    purifying_selection_factor = 0.2  # Ka is typically ~20% of Ks
    
    # The simple model's single prediction for Ka for any gene.
    simple_model_predicted_ka = baseline_mutation_rate * purifying_selection_factor
    
    print("--- Simulating the Challenge to Genetic Drift Models ---")
    print("\n[The Problem]: Models of genetic drift often use the synonymous substitution rate (Ks)")
    print("as a baseline 'neutral clock'. A major challenge arises when factors assumed to be neutral")
    print("and independent are actually correlated with selected factors.\n")

    # --- Step 2: Simulate the "Observed Reality" ---
    # In real genomes, local mutation rates vary. This means Ks is not constant,
    # and more importantly, it becomes correlated with Ka.
    np.random.seed(42) # For reproducibility
    # Let's simulate 5 genes where the local mutation rate varies.
    # This creates a range of observed Ks values.
    observed_ks_values = np.random.uniform(low=0.007, high=0.013, size=5)
    
    # In this correlated reality, the observed Ka depends on the observed Ks.
    # This is the key phenomenon.
    observed_ka_values = observed_ks_values * purifying_selection_factor + np.random.normal(0, 0.0001, 5)

    print("--- Comparing the Simple Model to Observed Reality ---")
    print(f"Simple Model's Universal Prediction for Ka: {simple_model_predicted_ka:.5f}")
    print("\nObserved Data (with correlation):")
    for i in range(len(observed_ks_values)):
        print(f"  Gene {i+1}: Observed Ks = {observed_ks_values[i]:.5f} -> Observed Ka = {observed_ka_values[i]:.5f}")

    # --- Step 3 & 4: Calculate Prediction Error to Demonstrate the Challenge ---
    # We will pick one gene to illustrate the model's failure.
    example_gene_index = 3 # Let's use the 4th gene
    example_observed_ka = observed_ka_values[example_gene_index]
    
    prediction_error = abs(simple_model_predicted_ka - example_observed_ka)
    
    print("\n--- Analysis of Prediction Failure ---")
    print("The simple model fails because it assumes a constant rate, but the reality is correlated.")
    print(f"For Gene {example_gene_index + 1}, the model predicts Ka = {simple_model_predicted_ka:.5f}, but the observed value is {example_observed_ka:.5f}.")
    print("\nThis discrepancy arises directly from the correlation between synonymous and nonsynonymous rates,")
    print("which is not accounted for in simple drift models. This directly supports Choice B.")

    # --- Step 5: Output the Final Equation with Numbers ---
    print("\nThe resulting error calculation is:")
    print("Final Equation: Prediction_Error = | Model_Prediction - Observed_Value |")
    print(f"Calculation for Gene {example_gene_index + 1}: Prediction_Error = | {simple_model_predicted_ka:.5f} - {example_observed_ka:.5f} | = {prediction_error:.5f}")

if __name__ == '__main__':
    genome_model_challenge()