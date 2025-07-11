import numpy as np

def grad_f(x):
    """
    Computes the gradient of f(x,y) = -x + 0.1y^2.
    x is a numpy array [x_coord, y_coord].
    """
    return np.array([-1.0, 0.2 * x[1]])

def project_to_tangent_cone(g, x_star):
    """
    Projects vector g onto the tangent cone of C at x_star.
    C = {(x,y) | x >= 0}.
    x_star = (0,0).
    Tangent cone T_C(x_star) is the right half-plane {(u,v) | u >= 0}.
    """
    # For x_star = (0,0), the tangent cone is the right half-plane.
    # Projection onto this set means the first component must be non-negative.
    u, v = g
    if u < 0:
        u = 0
    return np.array([u, v])

def is_in_normal_cone(v, x_star):
    """
    Checks if vector v is in the normal cone of C at x_star.
    C = {(x,y) | x >= 0}.
    x_star = (0,0).
    Normal cone N_C(x_star) is the non-positive x-axis {(u,v) | u <= 0, v = 0}.
    """
    # For x_star = (0,0), the normal cone is the negative x-axis including origin.
    u, v = v
    return u <= 0 and np.isclose(v, 0)

# Point to analyze
x_star = np.array([0.0, 0.0])

# --- Check for stationarity ---
# A point is stationary if Proj_{T_{x*}C}(-\nabla f(x*)) = 0
neg_grad_at_x_star = -grad_f(x_star)
v_star = project_to_tangent_cone(neg_grad_at_x_star, x_star)

is_stationary = np.allclose(v_star, 0)

print(f"Analyzing point x* = {x_star} for function f(x,y) = -x + 0.1y^2 on C = {{ (x,y) | x>=0 }}")
print(f"Negative gradient -\\nabla f(x*) = {-grad_f(x_star)}")
print(f"Projection onto tangent cone v* = Proj_T(- grad f(x*)) = {v_star}")
print(f"Is the point stationary? (v* == 0): {is_stationary}")
if not is_stationary:
    print("The point is NOT stationary.")

# --- Check if it is a fixed point of the iteration map ---
# A point is a fixed point if -v* is in the normal cone N_{x*}C
is_fixed_point = is_in_normal_cone(-v_star, x_star)

print(f"\nChecking if x* is a fixed point for Algorithm (2)...")
print(f"Vector to check for normal cone inclusion: -v* = {-v_star}")
print(f"Is -v* in the normal cone at x*? {is_fixed_point}")
if is_fixed_point:
    print("The point IS a fixed point of the iteration map.")

if not is_stationary and is_fixed_point:
    print("\nConclusion: The point is a non-stationary fixed point for Algorithm (2).")
    print("This demonstrates that it is possible for the algorithm to converge to a non-stationary point.")
