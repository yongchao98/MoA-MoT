import numpy as np

def run_dpgd_demonstration():
    """
    This function demonstrates that Doubly-projected Gradient Descent (DPGD),
    which is algorithm (2), can converge to a non-stationary point.

    We use the example discussed in the analysis:
    - Function: f(x) = x, defined on R^1.
    - Constraint set: C = R^1.
    - Gradient: grad(f)(x) = 1 everywhere.
    - Stationary points: None exist, since the gradient is never zero.
    
    The DPGD update simplifies to x_{k+1} = x_k - gamma_k.
    We choose a summable step-size sequence gamma_k = 0.5^k, which will cause
    the algorithm to converge to a limit point, which must be non-stationary.
    """

    print("--- Demonstration for Algorithm (2) ---")
    
    # Initial point
    x = 0.0

    # Number of iterations
    iterations = 50

    print(f"Running DPGD for f(x)=x with initial point x_0 = {x:.1f} and gamma_k = 0.5^k.")

    for k in range(iterations):
        gamma_k = 0.5**k
        # The gradient is always 1
        grad = 1.0
        x = x - gamma_k * grad

    # The analytical limit is x* = x_0 - sum_{k=0 to inf} 0.5^k = 0 - 2 = -2.
    x_star = -2.0
    
    print(f"\nAfter {iterations} iterations, the algorithm converges to x = {x:.8f}")
    print(f"The true analytical limit is x* = {x_star:.1f}")

    # Check the first-order stationarity condition at the limit point x*
    grad_at_limit = 1.0
    
    # The stationarity condition is ||Proj_{T_{x*}C}(-grad f(x*))|| = 0.
    # In this simplified case, it becomes ||-grad f(x*)|| = 0.
    norm_of_projected_gradient = np.linalg.norm(-grad_at_limit)
    
    print("\nChecking the stationarity condition at the limit point x* = -2.0:")
    
    print(f"The equation for stationarity is ||Proj_{{T_{{x*}}C}}(-\\nabla f(x*))|| = 0")
    print(f"Which for this example becomes || -({grad_at_limit}) || = 0")
    print("The numbers in the final equation are:")
    print(f"  Gradient at limit point: {grad_at_limit}")
    print(f"  Norm of projected gradient: {norm_of_projected_gradient}")
    
    print(f"\nThe equation is || -1.0 || = 1.0, which is not equal to 0.")
    print("\nConclusion: The algorithm has converged to a point that is NOT first-order stationary.")

run_dpgd_demonstration()