import numpy as np
from scipy.stats import pearsonr

def analyze_causation():
    """
    Analyzes the structural equation model E->A->B->C<-D<-E
    to determine if correlation between A and D implies causation.
    """
    print("Analyzing the Structural Equation Model: E->A->B->C<-D<-E")
    print("----------------------------------------------------------")

    # Step 1: Explain the causal structure regarding A and D.
    print("\nStep 1: Deconstructing the causal graph for A and D.")
    print("The model specifies the following causal paths originating from E:")
    print("1. E -> A (E causes A)")
    print("2. E -> D (E causes D)")
    print("This structure means that E is a 'common cause' for A and D.")
    print("There is no direct causal arrow pointing from A to D or from D to A.")

    # Step 2: Simulate the system using structural equations.
    print("\nStep 2: Simulating the system to demonstrate the relationship.")
    print("We can represent the causal effects with the following equations:")
    
    # Define and print the coefficients (the "numbers in the equation")
    coef_ea = 0.9
    coef_ed = 0.8
    print(f"Equation for A: A = {coef_ea} * E + (random noise_A)")
    print(f"Equation for D: D = {coef_ed} * E + (random noise_D)")
    print("Note: The other variables (B, C) are not needed to analyze the A-D relationship.")
    
    # Generate data
    n_samples = 20000
    E = np.random.normal(loc=0, scale=1, size=n_samples)
    noise_A = np.random.normal(loc=0, scale=0.5, size=n_samples)
    noise_D = np.random.normal(loc=0, scale=0.5, size=n_samples)
    
    A = coef_ea * E + noise_A
    D = coef_ed * E + noise_D

    # Step 3: Calculate the correlation from the simulated data.
    correlation, _ = pearsonr(A, D)
    print("\nStep 3: Calculating the correlation between A and D.")
    print(f"Based on {n_samples} simulated data points, the correlation is: {correlation:.4f}")

    # Step 4: Draw the final conclusion.
    print("\nStep 4: Conclusion.")
    print("The simulation confirms a strong, positive correlation between A and D.")
    print("However, this correlation arises because both A and D are caused by the confounding variable E.")
    print("An increase in E tends to increase both A and D, making them appear correlated.")
    print("This is a classic example of a spurious (non-causal) correlation due to confounding.")
    
    print("\n----------------------------------------------------------")
    print("Does correlation imply causation in this system? The final answer is:")
    print("No")

if __name__ == '__main__':
    analyze_causation()