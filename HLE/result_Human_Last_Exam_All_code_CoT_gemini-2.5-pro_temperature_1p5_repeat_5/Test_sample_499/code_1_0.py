import numpy as np

def illustrate_perturbation_scale():
    """
    Illustrates the importance of initialization magnitude in perturbation theory.

    Perturbation theory assumes the change in weights (dW) is small
    relative to the initial weights (W_initial). This code shows that
    for a fixed dW, the relative change is much larger when W_initial is small.
    """
    # Set a random seed for reproducibility
    np.random.seed(0)

    # A fixed perturbation (change in weights)
    dW = np.random.randn(5, 5) * 0.1
    norm_dW = np.linalg.norm(dW)

    print("--- Perturbation Details ---")
    print(f"Norm of the weight perturbation (||dW||): {norm_dW:.4f}\n")


    # Case 1: Small magnitude initialization
    # Scale factor for "small" initialization
    small_scale = 0.01
    W_initial_small = np.random.randn(5, 5) * small_scale
    norm_W_initial_small = np.linalg.norm(W_initial_small)

    # Calculate relative change
    relative_change_small = norm_dW / norm_W_initial_small

    print("--- Case 1: Small Initialization Magnitude ---")
    print(f"Initial Weight Norm (||W_initial_small||): {norm_W_initial_small:.4f}")
    print("Relative Change Equation: ||dW|| / ||W_initial_small||")
    print(f"Numbers for Equation: {norm_dW:.4f} / {norm_W_initial_small:.4f}")
    print(f"Resulting Relative Change: {relative_change_small:.4f}\n")
    print("This large relative change means a low-order perturbation theory is likely a poor approximation.\n")


    # Case 2: Large magnitude initialization
    # Scale factor for "large" initialization
    large_scale = 10.0
    W_initial_large = np.random.randn(5, 5) * large_scale
    norm_W_initial_large = np.linalg.norm(W_initial_large)

    # Calculate relative change
    relative_change_large = norm_dW / norm_W_initial_large

    print("--- Case 2: Large Initialization Magnitude ---")
    print(f"Initial Weight Norm (||W_initial_large||): {norm_W_initial_large:.4f}")
    print("Relative Change Equation: ||dW|| / ||W_initial_large||")
    print(f"Numbers for Equation: {norm_dW:.4f} / {norm_W_initial_large:.4f}")
    print(f"Resulting Relative Change: {relative_change_large:.4f}\n")
    print("This small relative change means a low-order perturbation theory is likely a good approximation (the 'lazy training' regime).")


if __name__ == '__main__':
    illustrate_perturbation_scale()
<<<D>>>