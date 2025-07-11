import math

def f(x):
    """A simple linear function."""
    return x

def grad_f(x):
    """Gradient of f(x)."""
    return 1.0

def proj_C(x):
    """Projection onto the constraint set C = [0, infinity)."""
    return max(0.0, x)

def proj_tangent_cone(x, v):
    """
    Projection onto the tangent cone of C at x.
    If x > 0, the tangent cone is R, so projection is identity.
    If x = 0, the tangent cone is [0, infinity).
    """
    if x > 0:
        return v
    else: # x == 0
        if v >= 0:
            return v
        else:
            return 0.0

def doubly_projected_gradient_descent_step(x_k, gamma_k):
    """Performs one step of Algorithm (2)."""
    # For our f(x)=x and x>0, T_x C = R, so Proj_T is identity.
    # The update simplifies to standard Projected Gradient Descent for x > 0.
    # x_{k+1} = Proj_C(x_k - gamma_k * grad_f(x_k))
    grad = grad_f(x_k)
    inner_proj = proj_tangent_cone(x_k, -grad)
    x_k_plus_1 = proj_C(x_k + gamma_k * inner_proj)
    return x_k_plus_1

def demonstrate_convergence_to_non_stationary():
    """
    Demonstrates that Algorithm (2) can converge to a non-stationary point
    by using a summable step-size sequence.
    """
    x_k = 5.0  # Initial point
    
    # For f(x)=x on C=[0,inf), the only stationary point is x=0.
    # Any other point is non-stationary.
    print(f"Starting at x_0 = {x_k}, which is not a stationary point.")

    # Use a summable step-size sequence
    total_gamma = 0
    num_iterations = 10000
    for k in range(1, num_iterations + 1):
        gamma_k = 1.0 / (k**2)
        total_gamma += gamma_k
        x_k = doubly_projected_gradient_descent_step(x_k, gamma_k)
    
    final_x = x_k
    
    # The theoretical limit is x_0 - sum(gamma_k)
    theoretical_limit = 5.0 - (math.pi**2 / 6)
    
    print(f"\nAfter {num_iterations} iterations with summable step sizes gamma_k = 1/k^2:")
    print(f"The algorithm converges to x = {final_x:.6f}")
    print(f"The theoretical limit is x_0 - pi^2/6 = {theoretical_limit:.6f}")
    
    # Check if the limit point is stationary
    grad_at_limit = grad_f(final_x)
    proj_grad_at_limit = proj_tangent_cone(final_x, -grad_at_limit)
    
    print(f"\nAt the limit point x = {final_x:.6f}:")
    print(f"The stationarity condition is ||Proj_T_x_C(-grad(f(x)))|| = 0.")
    print(f"-grad(f(x)) = {-grad_at_limit}")
    print(f"Proj_T_x_C(-grad(f(x))) = {proj_grad_at_limit}")
    
    is_stationary = (proj_grad_at_limit == 0)
    
    if not is_stationary:
        print("\nThe limit point is NOT a stationary point.")
    else:
        print("\nThe limit point IS a stationary point.")

demonstrate_convergence_to_non_stationary()