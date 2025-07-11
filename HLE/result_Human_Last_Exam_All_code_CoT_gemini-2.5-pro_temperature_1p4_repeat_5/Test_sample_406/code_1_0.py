import math

# This script illustrates the conditions for the existence of an FGH-tripled fixed point.
# It defines a set of functions F, G, H that satisfy the contractive conditions
# and then finds the fixed point using iterative approximation.

# --- Step 1: Define the functions F, G, H and their parameters ---
# The functions are defined as linear maps to easily control the Lipschitz constants.
# We will use the following function forms:
# F(x, y, z) = c_f + a1*x + a2*y + a3*z
# G(y, x, u) = c_g + b1*y + b2*x + b3*u  => G(y,x,y) = c_g + b2*x + (b1+b3)*y
# H(z, y, x) = c_h + g1*z + g2*y + g3*x

# Define constants (the constant part of the linear functions)
c_f, c_g, c_h = 1.0, 2.0, 3.0

# Define Lipschitz constants (the coefficients of the variables)
# For F(x,y,z):
a1, a2, a3 = 0.1, 0.1, 0.1
# For G(y,x,u):
b1, b2, b3 = 0.1, 0.1, 0.1
# For H(z,y,x):
g1, g2, g3 = 0.1, 0.1, 0.1

def F(x, y, z):
    """Function F: R*R*R -> R"""
    return c_f + a1 * x + a2 * y + a3 * z

def G(y, x, u):
    """Function G: R*R*R -> R"""
    return c_g + b1 * y + b2 * x + b3 * u

def H(z, y, x):
    """Function H: R*R*R -> R"""
    return c_h + g1 * z + g2 * y + g3 * x

# --- Step 2: Check if the contractive conditions are met ---
# For the operator T(x,y,z) = (F(x,y,z), G(y,x,y), H(z,y,x)) to be a contraction,
# the sum of Lipschitz constants for each variable dimension must be less than 1.
# (Here, X=Y=Z=R, and the metric is d(a,b)=|a-b|. The Lipschitz constants are the
# absolute values of the coefficients).

# Effective Lipschitz constants for the combined operator T
kx = abs(a1) + abs(b2) + abs(g3)
ky = abs(a2) + abs(b1) + abs(b3) + abs(g2)
kz = abs(a3) + abs(g1)

print("--- Checking Conditions for FGH-Tripled Fixed Point ---")
print("We consider spaces X, Y, Z to be the set of real numbers R (a complete metric space).")
print("The functions F, G, H are defined as linear maps with specific coefficients.")
print(f"Effective contraction factor for x-dimension (a1+b2+g3): {kx:.2f}")
print(f"Effective contraction factor for y-dimension (a2+b1+b3+g2): {ky:.2f}")
print(f"Effective contraction factor for z-dimension (a3+g1): {kz:.2f}")

k = max(kx, ky, kz)
print(f"\nThe overall contraction factor k = max({kx:.2f}, {ky:.2f}, {kz:.2f}) = {k:.2f}")

if k >= 1:
    print("\nThe conditions are not met, convergence is not guaranteed.")
else:
    print(f"\nSince k = {k:.2f} < 1, the conditions are met. A unique tripled fixed point exists.")
    print("-" * 50)

    # --- Step 3: Find the fixed point by iteration ---
    print("\n--- Finding the Fixed Point via Iteration ---")
    # Initial guess
    x_n, y_n, z_n = 0.0, 0.0, 0.0
    
    # Iteration parameters
    max_iterations = 100
    tolerance = 1e-9

    for i in range(max_iterations):
        # Apply the combined operator T(x,y,z) to get the next iterate
        x_n_plus_1 = F(x_n, y_n, z_n)
        # For G(y,x,y), the third argument 'u' is the same as the first argument 'y'
        y_n_plus_1 = G(y_n, x_n, y_n)
        z_n_plus_1 = H(z_n, y_n, x_n)
        
        # Check for convergence by measuring the distance between successive points
        error = math.sqrt((x_n_plus_1 - x_n)**2 + (y_n_plus_1 - y_n)**2 + (z_n_plus_1 - z_n)**2)
        
        # Update the point for the next iteration
        x_n, y_n, z_n = x_n_plus_1, y_n_plus_1, z_n_plus_1
        
        if error < tolerance:
            print(f"Convergence reached after {i+1} iterations.")
            break
    else:
        print("Maximum iterations reached without convergence.")
    
    # The converged point is our fixed point
    x_star, y_star, z_star = x_n, y_n, z_n

    print(f"\nThe calculated tripled fixed point is (x*, y*, z*) = ({x_star:.6f}, {y_star:.6f}, {z_star:.6f})")
    print("-" * 50)

    # --- Step 4: Verify the fixed point equations ---
    print("\n--- Verifying the Fixed Point Equations ---")
    
    # Calculate F(x*, y*, z*) and check if it equals x*
    f_val = F(x_star, y_star, z_star)
    print(f"Equation 1: F(x, y, z) = x")
    print(f"F({x_star:.6f}, {y_star:.6f}, {z_star:.6f}) = {f_val:.6f}")
    print(f"This should be equal to x* = {x_star:.6f}. Difference: {abs(f_val - x_star):.2e}")

    # Calculate G(y*, x*, y*) and check if it equals y*
    g_val = G(y_star, x_star, y_star)
    print(f"\nEquation 2: G(y, x, y) = y")
    print(f"G({y_star:.6f}, {x_star:.6f}, {y_star:.6f}) = {g_val:.6f}")
    print(f"This should be equal to y* = {y_star:.6f}. Difference: {abs(g_val - y_star):.2e}")

    # Calculate H(z*, y*, x*) and check if it equals z*
    h_val = H(z_star, y_star, x_star)
    print(f"\nEquation 3: H(z, y, x) = z")
    print(f"H({z_star:.6f}, {y_star:.6f}, {x_star:.6f}) = {h_val:.6f}")
    print(f"This should be equal to z* = {z_star:.6f}. Difference: {abs(h_val - z_star):.2e}")