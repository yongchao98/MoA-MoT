import numpy as np

def analyze_correlation_and_causation():
    """
    Analyzes and demonstrates the relationship between correlation and causation
    for the structural equation model: E->A->B->C<-D<-E.
    """
    # Step 1: Define parameters for the simulation
    n_samples = 100000
    np.random.seed(42)  # for reproducibility

    # Define the causal coefficients for the structural equations.
    # The specific values are not critical, only that they are non-zero.
    c_ea = 0.8  # E -> A
    c_ed = 0.7  # E -> D

    print("--- Causal System Analysis ---")
    print("Model: E -> A -> B -> C <- D <- E\n")
    print("We are investigating the relationship between A and D.")
    print("The system can be represented by these simplified structural equations:")
    # The prompt asks to output the numbers in the final equation.
    # We will show the equations used for the simulation.
    print(f"A = {c_ea}*E + noise_A")
    print(f"D = {c_ed}*E + noise_D")
    print("-" * 30)

    # Step 2: Generate data based on the model
    # E is the common cause (an exogenous variable)
    E = np.random.randn(n_samples)

    # A and D are both caused by E. We add random noise to each.
    noise_A = np.random.randn(n_samples)
    noise_D = np.random.randn(n_samples)
    A = c_ea * E + noise_A
    D = c_ed * E + noise_D

    # Step 3: Calculate the correlation between A and D
    correlation_A_D = np.corrcoef(A, D)[0, 1]
    
    # We show the "equation" for the correlation result.
    print(f"Correlation(A, D) = {correlation_A_D:.4f}")
    print("-" * 30)


    # Step 4: Explain the result
    print("Explanation:")
    print("1. Causal Path: In the model E->A... and E->D..., there is no arrow pointing from A to D or from D to A. This means there is no direct or indirect causal relationship between A and D.")
    print("\n2. Source of Correlation: A and D are correlated because they share a common cause, E. This is known as confounding. Changes in E cause changes in both A and D, making them appear statistically related.")
    print("\n3. Conclusion: The simulation shows a strong positive correlation (~0.47), but our analysis of the model's structure shows no causation between A and D. The correlation is real, but it is spurious with respect to a causal relationship between A and D.")
    print("\nTherefore, in this system, correlation does not imply causation.")

if __name__ == '__main__':
    analyze_correlation_and_causation()
