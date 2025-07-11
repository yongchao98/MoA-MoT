import numpy as np

def analyze_correlation_causation():
    """
    Simulates a structural equation model to demonstrate that correlation
    does not imply causation in the presence of a confounder.
    """
    # Set a seed for reproducibility of the random data
    np.random.seed(42)

    # Number of data points to simulate
    n_samples = 10000

    # Define the causal coefficients (the strength of each arrow in the model)
    c_ea = 1.0  # E -> A
    c_ed = 1.0  # E -> D
    c_ab = 0.8  # A -> B
    c_bc = 0.7  # B -> C
    c_dc = 0.6  # D -> C

    # --- Data Simulation ---
    # The model is E->A->B->C<-D<-E

    # E is an exogenous variable (no causes in the model), modeled with a normal distribution
    E = np.random.randn(n_samples)

    # A and D are both caused by the common variable E (the confounder)
    A = c_ea * E + np.random.randn(n_samples)
    D = c_ed * E + np.random.randn(n_samples)

    # B is caused by A
    B = c_ab * A + np.random.randn(n_samples)

    # C is caused by B and D (C is a "collider")
    C = c_bc * B + c_dc * D + np.random.randn(n_samples)

    # --- Analysis ---
    # Calculate the Pearson correlation coefficient between A and D
    correlation_matrix = np.corrcoef(A, D)
    correlation_ad = correlation_matrix[0, 1]

    # --- Output Results ---
    print("Causal Model: E -> A -> B -> C <- D <- E")
    print("\nAn observational study finds A and D are correlated.")
    print("The question is: Does this correlation imply causation between A and D?\n")

    print("We can simulate this system with the following structural equations:")
    print(f"A = {c_ea:.1f} * E + noise")
    print(f"D = {c_ed:.1f} * E + noise")
    # Note: A does not appear in the equation for D, and D does not appear in the equation for A.
    # This means there is no direct causation between them in our model.

    print(f"\nAfter simulating {n_samples} data points, we find:")
    print(f"The correlation between A and D is {correlation_ad:.4f}")

    print("\nConclusion:")
    print("A and D are indeed highly correlated. However, this correlation is not due to A causing D or D causing A.")
    print("It is because both variables are effects of the common cause E.")
    print("This is a classic example of confounding, where a third variable creates a spurious association between two others.")
    print("\nTherefore, in this system, correlation does NOT imply causation.")


if __name__ == "__main__":
    analyze_correlation_causation()

<<<No>>>