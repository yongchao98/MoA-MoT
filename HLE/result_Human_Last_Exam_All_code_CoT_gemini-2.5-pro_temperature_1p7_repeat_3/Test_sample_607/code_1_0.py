import numpy as np

def f(p):
    """Objective function f(x, y) = x + y."""
    return p[0] + p[1]

def grad_f(p):
    """Gradient of f(x, y)."""
    return np.array([1.0, 1.0])

def proj_C(p):
    """Projection onto C = union of x and y axes."""
    # Project onto x-axis
    p_x = np.array([p[0], 0])
    dist_x = np.linalg.norm(p - p_x)
    # Project onto y-axis
    p_y = np.array([0, p[1]])
    dist_y = np.linalg.norm(p - p_y)
    # Return the closer projection
    if dist_x <= dist_y:
        return p_x
    else:
        return p_y

def proj_tangent_cone(p, v):
    """Projection of vector v onto the tangent cone of C at point p."""
    # If p is not the origin, the tangent cone is the axis it lies on.
    # Note: A small tolerance is used for floating point comparison.
    if abs(p[1]) < 1e-9 and abs(p[0]) > 1e-9: # On x-axis, not origin
        return np.array([v[0], 0])
    elif abs(p[0]) < 1e-9 and abs(p[1]) > 1e-9: # On y-axis, not origin
        return np.array([0, v[1]])
    else: # At the origin (or very close), tangent cone is C itself.
        return proj_C(v)

def run_doubly_projected_gd():
    """Demonstrates convergence to a non-stationary point for algorithm (2)."""
    x = np.array([1.0, 0.0]) # Initial point
    gamma = 0.1 # Step size
    neg_grad = -grad_f(x) # Gradient is constant

    print(f"Running Doubly-Projected Gradient Descent...")
    print(f"f(x,y) = x+y, C = union of axes")
    print(f"Initial point x_0 = {x}")
    print(f"Algorithm converges to x* = (0,0)")
    
    for k in range(20):
        # The tangent cone projection step
        update_direction = proj_tangent_cone(x, neg_grad)
        # The update step
        x_plus_update = x + gamma * update_direction
        # The projection onto C
        x = proj_C(x_plus_update)
        # print(f"x_{k+1} = {x}") # Uncomment to see the iterates

    # After many iterations, x will be very close to the origin.
    # Let's check stationarity at the limit point x* = (0,0).
    x_star = np.array([0.0, 0.0])
    neg_grad_at_star = -grad_f(x_star)
    
    proj_on_tangent_at_star = proj_tangent_cone(x_star, neg_grad_at_star)
    
    norm_proj = np.linalg.norm(proj_on_tangent_at_star)

    print(f"\nLimit point is x* = {x_star.tolist()}")
    print(f"Is x* stationary? Check ||Proj_T_x*C(-grad f(x*))|| = 0")
    print(f"-grad f(x*) = {-neg_grad_at_star.tolist()}")
    print(f"T_x*C is C itself (union of axes).")
    print(f"Proj_T_x*C(-grad f(x*)) = {proj_on_tangent_at_star.tolist()}")
    print(f"Norm of projection = {norm_proj}")

    if norm_proj == 0:
        print("Conclusion: The limit point IS stationary.")
    else:
        print("Conclusion: The limit point IS NOT stationary.")

run_doubly_projected_gd()