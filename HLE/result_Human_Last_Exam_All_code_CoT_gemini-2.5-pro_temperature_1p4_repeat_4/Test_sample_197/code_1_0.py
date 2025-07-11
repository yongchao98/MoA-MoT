import numpy as np

def demonstrate_correlation_without_causation():
    """
    Simulates the structural equation model E->A->B->C<-D<-E to demonstrate
    that correlation between A and D does not imply causation.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # Number of data points to simulate
    n_samples = 10000

    # Define the structural equations with arbitrary coefficients.
    # We assume linear relationships for simplicity.
    # E is the common cause
    e = np.random.randn(n_samples)

    # A and D are caused by E
    # E -> A
    a = 0.8 * e + np.random.randn(n_samples) * 0.5
    # E -> D
    d = 0.7 * e + np.random.randn(n_samples) * 0.5

    # The rest of the model (not necessary for the A-D relationship, but included for completeness)
    # A -> B
    b = 0.6 * a + np.random.randn(n_samples) * 0.5
    # B -> C and D -> C (C is a collider)
    c = 0.5 * b + 0.4 * d + np.random.randn(n_samples) * 0.5

    # Calculate the correlation between A and D
    correlation = np.corrcoef(a, d)[0, 1]

    # Print the explanation and the result
    print("Analyzing the Structural Equation Model: E->A->B->C<-D<-E")
    print("-" * 60)
    print("1. Causal Structure: There is no direct or indirect causal path from A to D or D to A.")
    print("2. Common Cause: Both A and D are caused by the common variable E (a confounder).")
    print("3. Implication: This confounding structure induces a statistical correlation between A and D.")
    print("-" * 60)
    print("To demonstrate this, we simulated the system and calculated the correlation.")
    # The user requested to output the numbers in the final equation.
    # Here, the "equation" is the correlation calculation.
    print(f"Correlation(A, D) = {correlation:.4f}")
    print("\nConclusion: The high correlation is caused by the confounder E, not by a causal link between A and D.")
    print("Therefore, in this system, correlation does not imply causation.")


if __name__ == "__main__":
    demonstrate_correlation_without_causation()
