import numpy as np

def demonstrate_perturbation_theory():
    """
    Demonstrates how initial weight magnitude affects optimal parameters
    in a second-order approximation of a neural network's loss.
    """
    print("----------------------------------------------------------------------")
    print("Demonstration: Optimal Parameters Depend on Initialization Magnitude")
    print("----------------------------------------------------------------------")
    print("We analyze a simple neural network for a regression task.")
    print("Model: f(x, w) = w₂ * tanh(w₁ * x)")
    print("Loss is approximated by a quadratic: L(w) ≈ L(w₀) + gᵀ(w-w₀) + ½(w-w₀)ᵀH(w-w₀)")
    print("The optimal weights 'w_star' in this approximation are: w* ≈ w₀ - H⁻¹g")
    print("We will calculate w* for two different initial weight magnitudes.\n")

    # Data for a simple regression task
    X = np.array([-2, -1, 1, 2])
    Y_true = np.array([-1, -0.5, 0.5, 1])

    def network_model(x, w):
        """The neural network function."""
        w1, w2 = w
        return w2 * np.tanh(w1 * x)

    def jacobian(x, w):
        """Computes the Jacobian of the network output w.r.t. weights w."""
        w1, w2 = w
        # Derivative w.r.t w1
        d_dw1 = w2 * x / (np.cosh(w1 * x)**2)
        # Derivative w.r.t w2
        d_dw2 = np.tanh(w1 * x)
        return np.array([d_dw1, d_dw2])

    def analyze_initialization(w0):
        """
        Calculates gradient, Hessian, and optimal weights for a given initialization.
        """
        # Initialize gradient vector and Hessian matrix
        g = np.zeros(2)
        # Using the Gauss-Newton approximation for the Hessian: H ≈ JᵀJ
        H = np.zeros((2, 2))
        
        # Accumulate gradient and Hessian over all data points
        for x, y_true in zip(X, Y_true):
            J = jacobian(x, w0)
            y_pred = network_model(x, w0)
            residual = y_pred - y_true
            
            # Gradient of squared error loss is J * residual
            g += J * residual
            # Hessian (Gauss-Newton) is sum of outer products JᵀJ
            H += np.outer(J, J)
            
        # To ensure H is invertible, add a small regularization term (damping)
        H += np.eye(2) * 1e-5
        H_inv = np.linalg.inv(H)
        
        # The Newton step gives the change in weights
        delta_w = -H_inv @ g
        
        # The optimal parameters in the quadratic approximation
        w_star = w0 + delta_w
        
        # Print all the numbers in the final equation: w_star = w0 - H_inv @ g
        print(f"Initial weights w₀:\n{w0}")
        print("\nCalculated Hessian (H) at w₀:")
        print(np.round(H, 4))
        print("\nCalculated Gradient (g) at w₀:")
        print(np.round(g, 4))
        print(f"\nCalculated optimal weights w* ≈ w₀ - H⁻¹g:")
        print(f"{np.round(w_star, 4).tolist()} = {w0.tolist()} - H_inv @ {np.round(g, 4).tolist()}")
        return w_star

    # Case 1: Small magnitude initialization
    print("--- CASE 1: Small Magnitude Initialization ---")
    w_small = np.array([0.1, 0.1])
    w_star_small = analyze_initialization(w_small)
    print("\n")

    # Case 2: Large magnitude initialization
    print("--- CASE 2: Large Magnitude Initialization ---")
    w_large = np.array([2.0, 0.5]) # Kept w2 small to avoid extreme gradients
    w_star_large = analyze_initialization(w_large)
    print("\n")

    print("----------------------------------------------------------------------")
    print("Conclusion:")
    print(f"The calculated optimal parameters are different: {np.round(w_star_small, 4)} vs {np.round(w_star_large, 4)}.")
    print("This demonstrates that the optimal parameters under a second-order perturbation analysis")
    print("are determined by the properties of the Hessian and gradient at initialization, which are")
    print("themselves a function of the initial weight magnitudes.")
    print("----------------------------------------------------------------------")


if __name__ == "__main__":
    demonstrate_perturbation_theory()
<<<D>>>