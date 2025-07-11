import numpy as np

def solve_emi_shielding_force():
    """
    Calculates the force per unit area on a conducting plane in an EMI shielding setup.
    
    This function first prints the derived symbolic equation for the force. It then
    substitutes example numerical values for the physical parameters, prints the
    equation with these values, and finally computes and prints the numerical result.
    """
    
    # Print the derived symbolic formula
    print("The derived formula for the force per unit y-z area on the x = d interface is:")
    print("f/area = (μ₀/2) * (K₀² * sin²(a*y)) / [cosh(a*d) + (μ₀/μ) * sinh(a*d)]² * î_x\n")
    
    # Define example numerical values for the parameters
    K0 = 10.0  # Current sheet amplitude in A/m
    a = 15.0   # Spatial frequency in rad/m
    y = np.pi / (2 * a) # y-coordinate in m, chosen to maximize sin(a*y)
    d = 0.01   # Air gap thickness in m
    mu0 = 4 * np.pi * 1e-7  # Permeability of free space in H/m
    # Example relative permeability of the magnetic material (e.g., ferrite)
    mu_r = 1000.0
    mu = mu_r * mu0  # Permeability of the magnetic material in H/m

    # Calculate the components of the equation
    numerator = mu0 / 2 * K0**2 * np.sin(a * y)**2
    cosh_ad = np.cosh(a * d)
    sinh_ad = np.sinh(a * d)
    denominator = (cosh_ad + (mu0 / mu) * sinh_ad)**2
    
    force_magnitude = numerator / denominator

    # Print the equation with numerical values
    print("Using the following example values:")
    print(f"K₀ (current amplitude) = {K0} A/m")
    print(f"a (spatial frequency) = {a} rad/m")
    print(f"y (position) = {y:.4f} m (chosen such that sin(a*y) = 1)")
    print(f"d (air gap) = {d} m")
    print(f"μ₀ (vacuum permeability) = {mu0:.2e} H/m")
    print(f"μ (material permeability) = {mu:.2e} H/m (relative permeability μ_r = {mu_r})")
    print("\nThe equation with substituted values is:")
    print(f"f/area = ({mu0:.2e} / 2) * ({K0**2:.1f} * sin²({a*y:.2f})) / [cosh({a*d:.2f}) + ({mu0/mu:.1e}) * sinh({a*d:.2f})]² * î_x")
    print(f"f/area = ({numerator:.3e}) / [{cosh_ad:.3f} + ({mu0/mu:.1e}) * {sinh_ad:.3f}]² * î_x")
    print(f"f/area = ({numerator:.3e}) / ({denominator:.3f}) * î_x\n")

    # Print the final result
    print("The calculated force per unit area is:")
    print(f"f/area = {force_magnitude:.6f} N/m² * î_x")

# Execute the function to solve the problem
solve_emi_shielding_force()