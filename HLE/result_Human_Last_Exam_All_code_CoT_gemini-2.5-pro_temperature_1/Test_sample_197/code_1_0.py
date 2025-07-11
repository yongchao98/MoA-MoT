import pandas as pd
import numpy as np

def solve():
    """
    Analyzes and simulates the causal model E->A->B->C<-D<-E to determine
    if correlation between A and D implies causation.
    """
    # Set a seed for reproducibility
    np.random.seed(42)
    n_samples = 20000

    print("Analyzing the structural equation model: E->A->B->C<-D<-E\n")
    print("We will simulate this system to check the correlation between A and D.")

    # Define the causal coefficients for the linear equations
    # These represent the "strength" of the causal effects.
    # For example, in A = c_ea*E + error, c_ea is the coefficient.
    c_ea = 0.8
    c_ed = -0.7
    c_ab = 0.6
    c_bc = 0.5
    c_dc = 0.9

    print("\n--- Simulation Equations (with coefficients) ---")
    # E is an exogenous variable (it has no causes in the model)
    # We simulate it from a standard normal distribution.
    E = np.random.randn(n_samples)
    print("E = N(0, 1)")

    # A and D are caused by E
    A = c_ea * E + np.random.randn(n_samples) * 0.5
    D = c_ed * E + np.random.randn(n_samples) * 0.5
    print(f"A = {c_ea}*E + error_A")
    print(f"D = {c_ed}*E + error_D")

    # B is caused by A
    B = c_ab * A + np.random.randn(n_samples) * 0.5
    print(f"B = {c_ab}*A + error_B")

    # C is caused by B and D (C is a collider)
    C = c_bc * B + c_dc * D + np.random.randn(n_samples) * 0.5
    print(f"C = {c_bc}*B + {c_dc}*D + error_C")
    print("------------------------------------------------\n")


    # Create a pandas DataFrame to hold the simulated data
    df = pd.DataFrame({'A': A, 'D': D})

    # Calculate the correlation between A and D
    correlation_ad = df['A'].corr(df['D'])

    print(f"The simulated correlation between A and D is: {correlation_ad:.4f}")

    print("\nExplanation:")
    print("The simulation shows a strong correlation between A and D.")
    print("However, inspecting the model reveals no causal path from A to D or D to A.")
    print("The correlation is entirely due to the common cause (confounder) E.")
    print("Therefore, in this system, correlation does not imply causation.")
    print("\nFinal Answer:")

solve()
print("<<<No>>>")