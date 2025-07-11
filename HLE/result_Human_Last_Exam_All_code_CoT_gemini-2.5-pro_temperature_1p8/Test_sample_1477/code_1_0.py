import numpy as np

def demonstrate_convergence_bound():
    """
    This function demonstrates the convergence condition for gradient descent
    on a 1D logistic regression problem. It shows that for a learning rate
    gamma < 2/L, the algorithm converges, while for gamma > 2/L, it can diverge.
    """
    
    # --- Function Definitions ---
    
    def sigma(t):
        """Sigmoid function."""
        # Clip to avoid overflow in exp
        t = np.clip(t, -500, 500)
        return 1 / (1 + np.exp(-t))

    def R_prime(w, x, y):
        """Computes the gradient R'(w)."""
        N = len(x)
        # Note: The formula from the derivation is R'(w) = (1/N) * sum(-sigma(-y*w*x) * y * x)
        # Let's use this correct version.
        grad = np.sum(-sigma(-y * w * x) * y * x) / N
        return grad

    def R_pprime(w, x, y):
        """Computes the second derivative R''(w)."""
        N = len(x)
        hess = np.sum(x**2 * sigma(y * w * x) * sigma(-y * w * x)) / N
        return hess

    # --- Problem Setup ---
    # Define a simple non-separable dataset in 1D
    x = np.array([1.0, 2.0])
    y = np.array([1.0, -1.0])
    N = len(x)

    print("--- Logistic Regression Gradient Descent Analysis ---")
    print(f"Using dataset: x = {x}, y = {y}\n")
    
    # --- Calculate L and M ---
    
    # L is the maximum of the second derivative, which occurs at w = 0
    L = R_pprime(0, x, y)
    
    # The theoretical largest upper bound for the learning rate is M = 2/L
    M = 2 / L

    print(f"1. Calculating the smoothness constant L:")
    L_equation_val = (1/(4*N)) * np.sum(x**2)
    print(f"   L = sup R''(w) = R''(0) = (1/(4*N)) * sum(x_i^2)")
    print(f"   L = (1 / (4 * {N})) * ({x[0]**2} + {x[1]**2}) = {L_equation_val:.4f}")
    
    print("\n2. Calculating the upper bound for the learning rate M:")
    print(f"   The convergence condition is gamma < 2/L.")
    print(f"   Therefore, M = 2 / L")
    print(f"   M = 2 / {L:.4f} = {M:.4f}\n")
    
    print("--- Running Gradient Descent Simulations ---")

    def gradient_descent(gamma, w_init, x, y, n_iter=20):
        """Performs gradient descent and shows the results."""
        w = w_init
        print(f"\nRunning with learning rate gamma = {gamma:.4f} (Starting at w={w_init:.2f})")
        print("-" * 50)
        print(f"{'Iteration':<10} | {'w_k':<20} | {'|w_k - w_{k-1}|':<20}")
        print("-" * 50)
        
        is_diverging = False
        for i in range(n_iter):
            w_prev = w
            grad = R_prime(w, x, y)
            w = w - gamma * grad
            
            # Print status
            if i < 10 or i == n_iter - 1:
                print(f"{i+1:<10} | {w:<20.6f} | {abs(w - w_prev):<20.6f}")

            # Check for divergence
            if np.isnan(w) or np.abs(w) > 1e6:
                print("\n   Divergence detected!")
                is_diverging = True
                break
        
        if not is_diverging:
             print("\n   Converged.")
        print("-" * 50)


    # Case 1: gamma < M (should converge)
    gamma_converge = M - 0.1
    # We start near w=0, where the curvature is highest, to test the bound robustly.
    w_init = 0.1 
    gradient_descent(gamma_converge, w_init, x, y)

    # Case 2: gamma > M (expected to diverge)
    gamma_diverge = M + 0.1
    gradient_descent(gamma_diverge, w_init, x, y)

# Execute the demonstration
demonstrate_convergence_bound()