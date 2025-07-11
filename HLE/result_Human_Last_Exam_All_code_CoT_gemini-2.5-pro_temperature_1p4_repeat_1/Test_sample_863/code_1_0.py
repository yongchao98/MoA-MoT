import sys

def solve_for_chi_star(chi, Nx):
    """
    Calculates chi* based on the given chi and the demagnetizing factor Nx.
    
    The problem is to find chi* that satisfies the equation:
    Nm(a/b, chi) + Nm(b/a, chi*) = 1
    
    where Nm is the magnetometric demagnetizing factor, given by:
    Nm = N / (1 + N * chi)
    
    For an infinitely long prism, the fluxmetric demagnetizing factors satisfy:
    Nx + Ny = 1
    
    The equation can be written as:
    Nx / (1 + Nx * chi) + Ny / (1 + Ny * chi*) = 1
    
    We solve this equation for chi*.
    Let's rearrange the equation:
    Ny / (1 + Ny * chi*) = 1 - Nx / (1 + Nx * chi)
    Ny / (1 + Ny * chi*) = (1 + Nx * chi - Nx) / (1 + Nx * chi)
    Since Ny = 1 - Nx, the right side becomes:
    Ny / (1 + Ny * chi*) = (Ny + Nx * chi) / (1 + Nx * chi)
    
    Now, solving for chi*:
    1 + Ny * chi* = Ny * (1 + Nx * chi) / (Ny + Nx * chi)
    Ny * chi* = (Ny * (1 + Nx * chi) / (Ny + Nx * chi)) - 1
    Ny * chi* = (Ny + Ny * Nx * chi - Ny - Nx * chi) / (Ny + Nx * chi)
    Ny * chi* = (Ny * Nx * chi - Nx * chi) / (Ny + Nx * chi)
    Ny * chi* = chi * Nx * (Ny - 1) / (Ny + Nx * chi)
    Since Ny - 1 = -Nx:
    Ny * chi* = chi * Nx * (-Nx) / (Ny + Nx * chi)
    chi* = - (Nx**2 * chi) / (Ny * (Ny + Nx * chi))
    
    This is the general formula. Let's implement it.
    """
    
    # Check for invalid inputs
    if not -1 < Nx < 1:
        print("Error: Nx must be between 0 and 1.", file=sys.stderr)
        return None
        
    Ny = 1.0 - Nx

    # Denominator of the expression for chi*
    # Handle the case where the denominator might be zero
    denominator = Ny * (Ny + Nx * chi)
    if denominator == 0:
        print("Error: Calculation resulted in division by zero.", file=sys.stderr)
        return None
        
    # Numerator of the expression for chi*
    numerator = - (Nx**2 * chi)
    
    chi_star = numerator / denominator
    
    # The final equation is chi* = numerator / denominator. Let's output all parts.
    print(f"Given values:")
    print(f"chi = {chi}")
    print(f"N_x = {Nx}")
    print(f"Derived N_y = 1 - N_x = {Ny}\n")
    print(f"Calculation of chi*:")
    print(f"Numerator = - (N_x^2 * chi) = {numerator}")
    print(f"Denominator = N_y * (N_y + N_x * chi) = {denominator}")
    print(f"Result: chi* = {chi_star}")
    
    return chi_star

# --- Example Usage ---
# The demagnetizing factor Nx depends on the aspect ratio a/b.
# For a square cross-section (a=b), Nx = Ny = 0.5.
# Let's use this special case as an example.
print("--- Example: Square Cross-Section (a=b) ---")
# Let's assume a material susceptibility
chi_example = 10.0
Nx_example = 0.5

# For this special case, the formula simplifies:
# chi* = - (0.5^2 * chi) / (0.5 * (0.5 + 0.5 * chi))
# chi* = - (0.25 * chi) / (0.25 * (1 + chi))
# chi* = -chi / (1 + chi)
# For chi=10, chi* should be -10 / 11 = -0.9090...
print(f"Using chi = {chi_example}, Nx = {Nx_example}")
solve_for_chi_star(chi_example, Nx_example)

# Another example: a rectangle where Nx is different from 0.5
print("\n--- Example: Rectangular Cross-Section ---")
chi_example_2 = 20.0
Nx_example_2 = 0.7 
print(f"Using chi = {chi_example_2}, Nx = {Nx_example_2}")
solve_for_chi_star(chi_example_2, Nx_example_2)
