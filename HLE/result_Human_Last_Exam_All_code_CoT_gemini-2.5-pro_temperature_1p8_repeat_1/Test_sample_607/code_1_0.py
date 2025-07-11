import numpy as np

def demonstrate_convergence_to_non_stationary_point():
    """
    This function demonstrates that Doubly-projected gradient descent (2)
    can converge to a non-stationary point given an appropriate choice
    of step sizes gamma_k.
    """
    print("Analyzing Algorithm (2): Doubly-projected gradient descent")
    print("------------------------------------------------------------")

    # Define the problem: f(x) = x on C = R.
    # Gradient is nabla f(x) = 1 for all x.
    # A stationary point x* must satisfy ||Proj_{T_{x*}C}(-nabla f(x*))|| = 0.
    # For C=R and f(x)=x, this becomes ||-1|| = 1 = 0, which is impossible.
    # Thus, no stationary points exist for this problem.
    # If the algorithm converges to a point, it must be non-stationary.
    
    # Let's set up the simulation for algorithm (2)
    # The update rule is x_{k+1} = x_k - gamma_k
    x0 = 10.0
    q = 0.5 # For gamma_k = q^k
    num_iterations = 50

    # The series sum(q^k) from k=0 to infinity converges to 1/(1-q)
    series_sum = 1 / (1 - q)
    theoretical_limit = x0 - series_sum

    print(f"Problem setup: f(x) = x, C = R, x_0 = {x0}")
    print(f"Step sizes: gamma_k = ({q})^k")
    print(f"The sequence x_k should converge to x* = x_0 - sum(gamma_k).")
    print(f"Theoretical Limit Equation: x* = {x0} - {series_sum} = {theoretical_limit}")
    print("------------------------------------------------------------")
    print("Running simulation...")

    x_k = x0
    print(f"x_0 = {x_k:.8f}")
    for k in range(num_iterations):
        gamma_k = q**k
        x_k_next = x_k - gamma_k
        if k < 5 or k > num_iterations - 5:
            print(f"k={k:2d}, gamma_k={gamma_k:.6f}, x_{k+1}={x_k_next:.8f}")
        x_k = x_k_next

    nabla_f_at_limit = 1.0
    stationarity_norm = np.linalg.norm(-nabla_f_at_limit)
    
    print("------------------------------------------------------------")
    print("Simulation finished.")
    print(f"Final iterate value x_{num_iterations} = {x_k:.8f}")
    print(f"This is close to the theoretical limit x* = {theoretical_limit:.8f}.")
    print("\nChecking the stationarity of the limit point x*:")
    print(f"Gradient at x* is nabla f(x*) = {nabla_f_at_limit}")
    print(f"Stationarity measure ||Proj_T(-nabla f(x*))|| = {stationarity_norm}")
    print("\nSince the stationarity measure is not 0, the limit point is NOT stationary.")
    print("This shows that algorithm (2) can converge to a non-stationary point.")

if __name__ == '__main__':
    demonstrate_convergence_to_non_stationary_point()