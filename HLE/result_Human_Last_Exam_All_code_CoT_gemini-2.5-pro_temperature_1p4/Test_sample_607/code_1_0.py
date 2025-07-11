import numpy as np

def run_dpgd_on_linear_function():
    """
    This function demonstrates that Doubly-Projected Gradient Descent (DPGD) can converge
    to a non-stationary point if a summable step-size sequence is used.

    We consider the simple 1D function f(x) = -x over the constraint set C = R.
    The gradient is f'(x) = -1 for all x.
    A stationary point must satisfy the condition grad(f(x)) = 0, which is -1 = 0,
    an impossibility. Therefore, no stationary points exist for this function.

    We will show that DPGD converges to a specific point, which must therefore be
    a non-stationary point.
    """
    # Initial point
    x = 0.0
    
    # Parameters
    num_iterations = 50
    # A summable step-size sequence, gamma_k = 0.5^(k+1)
    # The sum of this series is 1.0
    gammas = [0.5**(k + 1) for k in range(num_iterations)]

    # The function f(x) = -x has a constant gradient f'(x) = -1
    grad_f = -1.0
    
    # Since C = R, the tangent cone is R, and the projection onto it is the identity.
    # The outer projection onto C is also the identity.
    # The DPGD update simplifies to: x_{k+1} = x_k + gamma_k * (-grad_f)
    projected_neg_grad = -grad_f
    
    print(f"Running DPGD for f(x) = -x with x_0 = {x} and summable step sizes.")
    
    for k in range(num_iterations):
        x = x + gammas[k] * projected_neg_grad

    limit_point = x
    # The theoretical limit is x_0 + Sum(gammas) * (-grad_f) = 0 + 1.0 * 1.0 = 1.0
    
    print(f"The algorithm converged to the point x* = {limit_point:.6f}")
    
    # Check the stationarity condition at the limit point
    grad_at_limit = -1.0
    stationarity_condition_norm = np.linalg.norm(grad_at_limit)
    
    print(f"The gradient at this point is grad(f(x*)) = {grad_at_limit}")
    print(f"The norm of the projected negative gradient (stationarity condition) is ||Proj(-grad(f(x*)))|| = {stationarity_condition_norm:.6f}")
    
    if stationarity_condition_norm > 1e-6:
        print("Since the stationarity condition norm is not zero, the point is not first-order stationary.")
    else:
        print("The point is first-order stationary.")

run_dpgd_on_linear_function()