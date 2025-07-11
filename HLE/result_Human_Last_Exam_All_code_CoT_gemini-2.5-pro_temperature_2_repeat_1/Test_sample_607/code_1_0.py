import numpy as np

def f(x):
    """The function to minimize."""
    return x

def grad_f(x):
    """The gradient of the function."""
    return 1.0

def proj_C(x):
    """Projection onto the constraint set C = [0, infinity)."""
    return np.maximum(0, x)

def proj_TxC(v, x):
    """Projection onto the tangent cone of C at x."""
    if x > 0:
        # For x in the interior of C, the tangent cone is R^1.
        return v
    elif x == 0:
        # For x at the boundary, the tangent cone is [0, infinity).
        return np.maximum(0, v)
    else:
        # x is not in C, this case should not be reached with this algorithm.
        raise ValueError("x is not in the constraint set C")

def check_stationarity(x):
    """Check the first-order stationarity condition at point x."""
    neg_grad = -grad_f(x)
    projected_neg_grad = proj_TxC(neg_grad, x)
    # The condition is ||proj|| = 0, which is equivalent to proj = 0 in 1D.
    return abs(projected_neg_grad)

def algorithm_2(x0, iterations):
    """
    Implements the doubly-projected gradient descent algorithm.
    """
    x = x0
    print(f"Initial point x_0 = {x}")
    
    # Pathological step sizes with a finite sum
    gammas = [(0.5)**(k + 1) for k in range(iterations)]
    
    for k in range(iterations):
        gamma_k = gammas[k]
        
        # Calculate negative gradient
        neg_grad = -grad_f(x)
        
        # Project negative gradient onto tangent cone
        d = proj_TxC(neg_grad, x)
        
        # Update step
        x_next = proj_C(x + gamma_k * d)
        
        x = x_next
        if k < 5 or k == iterations -1:
            print(f"Iteration {k+1:4d}: x_{k+1} = {x:.8f}")

    print("\n--- Final Analysis ---")
    print(f"The algorithm converged to a point near: {x}")
    stationarity_metric = check_stationarity(x)
    print(f"Stationarity metric at the final point x = {x:.4f} is ||Proj_T_xC(-grad f(x))|| = {stationarity_metric:.4f}")

    if stationarity_metric > 1e-6:
        print("The final point is NOT a first-order stationary point.")
    else:
        print("The final point IS a first-order stationary point.")

# --- Run the demonstration ---
initial_point = 2.0
num_iterations = 30
algorithm_2(initial_point, num_iterations)
