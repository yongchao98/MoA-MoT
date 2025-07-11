import numpy as np

def analyze_causal_system():
    """
    Simulates the causal system E->A->B->C<-D<-E to check the correlation
    between A and D and explain the implications for causality.
    """
    # Set a seed for reproducibility of the random data
    np.random.seed(42)

    # Number of data points in our simulated observational study
    n_samples = 20000

    # Step 1: Define the Structural Equation Model (SEM)
    # We model the relationships using linear equations with random noise.
    # The coefficients (e.g., 0.7, 0.6) represent the strength of the causal effect.

    # E is a root cause, so we generate it from a standard normal distribution.
    E = np.random.randn(n_samples)

    # A is caused by E. The equation is: A = 0.7*E + noise
    A = 0.7 * E + np.random.randn(n_samples) * 0.5

    # D is caused by E. The equation is: D = 0.6*E + noise
    D = 0.6 * E + np.random.randn(n_samples) * 0.5

    # B is caused by A. The equation is: B = 0.8*A + noise
    B = 0.8 * A + np.random.randn(n_samples) * 0.5

    # C is caused by B and D. The equation is: C = 0.5*B + 0.4*D + noise
    C = 0.5 * B + 0.4 * D + np.random.randn(n_samples) * 0.5

    # Step 2: Calculate the correlation between A and D
    # In our model, there is no direct arrow A->D or D->A.
    # Any correlation found would be due to other factors.
    correlation_matrix = np.corrcoef(A, D)
    correlation_ad = correlation_matrix[0, 1]

    # Step 3: Print the results and the final conclusion
    print("Analysis of the Causal System: E->A->B->C<-D<-E")
    print("=" * 55)
    print(f"Finding: The calculated correlation between A and D is {correlation_ad:.4f}")
    print("\nExplanation:")
    print("1. In this model, there is NO direct or indirect causal path from A to D or from D to A.")
    print("2. The variables A and D share a common cause, E (since E->A and E->D).")
    print("3. This common cause, E, is a confounder that creates a statistical correlation between A and D.")
    print("4. The correlation is real, but it is not the result of a causal relationship between A and D.")
    print("\nFinal Answer Formulation:")
    print("The question is: Does correlation imply causation in this system?")
    print("The answer is based on the fact that the correlation between A and D is entirely explained by the confounder E.")

if __name__ == '__main__':
    analyze_causal_system()
