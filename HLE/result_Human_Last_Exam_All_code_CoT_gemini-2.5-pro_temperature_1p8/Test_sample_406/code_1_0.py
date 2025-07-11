# Step 1: Define example functions F, G, and H that form a contraction mapping.
# We'll use simple linear functions for this demonstration.
def F(x, y, z):
    """Function F: X*Y*Z -> X"""
    # The sum of absolute values of coefficients for x,y,z (0.1+0.1+0.1=0.3) is less than 1.
    return 0.1 * x + 0.1 * y - 0.1 * z + 5.0

def G(y, x, y_prime):
    """Function G: Y*X*Y -> Y. Note y_prime is the same as y in a fixed point context."""
    return 0.2 * y - 0.1 * x + 0.1 * y_prime + 2.0

def H(z, y, x):
    """Function H: Z*Y*X -> Z"""
    return 0.1 * z + 0.2 * y - 0.1 * x + 1.0

# Step 2: Initialize a starting point (x0, y0, z0)
x, y, z = 0.0, 0.0, 0.0

# Step 3: Iteratively apply the functions to find the fixed point.
# Because the mapping is a contraction, this process will converge to the unique fixed point.
num_iterations = 100
for i in range(num_iterations):
    x_next = F(x, y, z)
    # For G(y, x, y), we use the current 'y' for both y inputs.
    y_next = G(y, x, y)
    z_next = H(z, y, x)
    
    # Update the current point
    x, y, z = x_next, y_next, z_next

# Step 4: The final (x, y, z) is our approximate tripled fixed point.
# Print the results and verify them by substituting back into the equations.
print(f"Found the approximate tripled fixed point after {num_iterations} iterations:")
print(f"(x, y, z) = ({x:.6f}, {y:.6f}, {z:.6f})")
print("-" * 40)
print("Verifying the fixed point equations:")
print("-" * 40)

# Check F(x,y,z) = x
fx_val = F(x, y, z)
print(f"F({x:.6f}, {y:.6f}, {z:.6f}) = {fx_val:.6f}")
print(f"This should be equal to x = {x:.6f}\n")

# Check G(y,x,y) = y
gy_val = G(y, x, y)
print(f"G({y:.6f}, {x:.6f}, {y:.6f}) = {gy_val:.6f}")
print(f"This should be equal to y = {y:.6f}\n")

# Check H(z,y,x) = z
hz_val = H(z, y, x)
print(f"H({z:.6f}, {y:.6f}, {x:.6f}) = {hz_val:.6f}")
print(f"This should be equal to z = {z:.6f}")