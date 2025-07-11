import numpy as np

def f(x):
    """The objective function."""
    return 0.5 * x**2

def grad_f(x):
    """The gradient of the objective function."""
    return x

def run_doubly_projected_gd_counterexample():
    """
    This function demonstrates that Doubly-Projected Gradient Descent (2)
    can converge to a non-stationary point.

    We use f(x) = 0.5*x^2 and the constraint set C = [10, infinity).
    The only stationary point is x=10.

    We show that the algorithm can be made to converge to x=20.
    """
    C_lower_bound = 10.0
    stationary_point = 10.0
    target_point = 20.0
    
    print(f"Objective function f(x) = 0.5 * x^2")
    print(f"Constraint set C = [{C_lower_bound}, infinity)")
    print(f"The unique stationary point is x = {stationary_point}")
    print(f"We will show convergence to the non-stationary point x = {target_point}\n")

    # Start at k=1
    x_k = target_point + 1.0 / 1.0 
    
    print("Running simulation...")
    print("k |       x_k       |      gamma_k")
    print("------------------------------------------")
    
    # We print the state for k=1 manually
    gamma_k = 1.0 / (20.0 * 1**2 + 21.0 * 1 + 1.0)
    print(f"{1:^1} | {x_k:^15.10f} | {gamma_k:^15.10f}")

    # Iterate for k from 2 up to N
    num_iterations = 15
    for k in range(1, num_iterations + 1):
        # For our chosen sequence x_k = 20 + 1/k, all points are in the interior of C.
        # So, the tangent cone T_{x_k}C is the entire real line R.
        # The projection onto the tangent cone is the identity.
        # Proj_{T_{x_k}C}(-grad_f(x_k)) = -grad_f(x_k) = -x_k
        
        # We define gamma_k specifically to generate our desired sequence
        # x_k = 20 + 1/k. The update rule is x_{k+1} = (1-gamma_k)*x_k.
        # This gives gamma_k = 1 - x_{k+1}/x_k.
        # Substituting x_k=20+1/k and x_{k+1}=20+1/(k+1) gives:
        gamma_k = 1.0 / (20.0 * k**2 + 21.0 * k + 1.0)
        
        # Calculate the update
        update_direction = -grad_f(x_k) # Projection on T_xk C is identity
        y_k = x_k + gamma_k * update_direction
        
        # Project back onto the constraint set C
        x_k_plus_1 = np.maximum(C_lower_bound, y_k)
        
        # For the final report
        if k < num_iterations:
            print(f"{k+1:^1} | {x_k_plus_1:^15.10f} | {1.0 / (20.0 * (k+1)**2 + 21.0 * (k+1) + 1.0):^15.10f}")
        
        # Update x_k for the next iteration
        x_k = x_k_plus_1

    print("------------------------------------------")
    print(f"\nAfter {num_iterations} iterations, the algorithm is near x = {x_k:.4f}.")
    print(f"This is close to our target of {target_point}, not the stationary point {stationary_point}.")
    print("This demonstrates that algorithm (2) can converge to a non-stationary point.")
    
    # We can also explicitly print the final equation.
    # The final point x_k should be approximately equal to the target point 20.
    # The stationary condition at x=20 is ||Proj_{T_{20}C}(-grad(f(20)))|| = ||-20|| = 20 != 0
    grad_at_limit = grad_f(target_point)
    proj_norm_at_limit = np.abs(-grad_at_limit)
    print(f"\nFinal check: At the limit point x* = {target_point}:")
    print(f"The stationarity condition is ||Proj(T_x* C, -grad(f(x*)))|| = 0.")
    print(f"For x* = {target_point}, this is ||Proj(R, -{int(target_point)})|| = ||-{int(target_point)}|| = {int(proj_norm_at_limit)}.")
    print(f"Since {int(proj_norm_at_limit)} is not equal to 0, the point is not stationary.")


if __name__ == '__main__':
    run_doubly_projected_gd_counterexample()