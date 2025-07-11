import math

def evaluate_1s_nuclear_attraction(zeta_i, zeta_j):
    """
    Calculates the nuclear attraction integral <phi_i | 1/r | phi_j>
    for two 1s Slater-Type Orbitals centered at the origin.

    Args:
        zeta_i (float): The orbital exponent of the first 1s STO.
        zeta_j (float): The orbital exponent of the second 1s STO.

    Returns:
        float: The value of the integral.
    """
    if zeta_i <= 0 or zeta_j <= 0:
        raise ValueError("Orbital exponents (zeta) must be positive.")
    
    numerator = 4 * (zeta_i * zeta_j)**1.5
    denominator = (zeta_i + zeta_j)**2
    return numerator / denominator

def explain_and_calculate():
    """
    Prints the derivation and example calculations.
    """
    print("### Evaluating the integral <phi_i | 1/r | phi_j> for 1s Slater Orbitals ###")
    
    print("\n--- Step 1: Define the Orbitals ---")
    print("A 1s Slater-Type Orbital (STO) centered at the origin is given by the function:")
    print("phi(r; zeta) = N * exp(-zeta*r), where N is the normalization constant.")
    print("The normalized form is: phi(r; zeta) = (zeta**3 / pi)**0.5 * exp(-zeta*r)")
    print("We consider two such orbitals, phi_i and phi_j, with exponents zeta_i and zeta_j, both centered at the origin.")
    
    print("\n--- Step 2: Set up the Integral ---")
    print("The integral to evaluate is: I_ij = Integral( phi_i^*(r) * (1/r) * phi_j(r) dV ) over all space.")
    print("Since the functions are real, phi_i^* = phi_i. Substituting the STO expressions:")
    print("I_ij = Integral( (zeta_i**3/pi)**0.5 * exp(-zeta_i*r) * (1/r) * (zeta_j**3/pi)**0.5 * exp(-zeta_j*r) dV )")
    print("Combining constant terms gives:")
    print("I_ij = ((zeta_i*zeta_j)**1.5 / pi) * Integral( (1/r) * exp(-(zeta_i + zeta_j)*r) dV )")
    
    print("\n--- Step 3: Integrate using Spherical Coordinates ---")
    print("In spherical coordinates, the volume element is dV = r**2 * sin(theta) dr dtheta dphi.")
    print("The integrand does not depend on the angles theta and phi. The integral over the solid angle is 4*pi.")
    print("I_ij = ((zeta_i*zeta_j)**1.5 / pi) * 4*pi * Integral from 0 to infinity of [ (1/r) * exp(-(zeta_i + zeta_j)*r) * r**2 dr ]")
    print("This simplifies to:")
    print("I_ij = 4 * (zeta_i*zeta_j)**1.5 * Integral from 0 to infinity of [ r * exp(-(zeta_i + zeta_j)*r) dr ]")
    
    print("\n--- Step 4: Solve the Radial Integral ---")
    print("The radial integral is a standard gamma function form: Integral( x^n * exp(-a*x) dx ) = n! / a**(n+1)")
    print("In our case, n=1, x=r, and a = (zeta_i + zeta_j).")
    print("So, the integral evaluates to: 1 / (zeta_i + zeta_j)**2")

    print("\n--- Step 5: The Final Analytical Formula ---")
    print("Substituting the result of the radial integral back into the expression for I_ij gives the final formula:")
    print("I_ij = 4 * (zeta_i * zeta_j)**1.5 / (zeta_i + zeta_j)**2")
    print("-" * 60)
    
    # --- Example Calculations ---
    print("\n### Example Calculations ###")
    
    # Example 1: Diagonal case (i=j), representing <1/r> for a hydrogen atom ground state.
    zeta_1_i = 1.0
    zeta_1_j = 1.0
    result_1 = evaluate_1s_nuclear_attraction(zeta_1_i, zeta_1_j)
    
    print(f"\n1. Diagonal Case (phi_i = phi_j): Let zeta_i = {zeta_1_i}, zeta_j = {zeta_1_j}")
    print(f"   I = 4 * ({zeta_1_i} * {zeta_1_j})**1.5 / ({zeta_1_i} + {zeta_1_j})**2")
    print(f"   I = 4 * ({zeta_1_i * zeta_1_j})**1.5 / ({zeta_1_i + zeta_1_j})**2")
    print(f"   I = {4 * (zeta_1_i * zeta_1_j)**1.5} / { (zeta_1_i + zeta_1_j)**2 } = {result_1}")
    
    # Example 2: Off-diagonal case (i != j)
    zeta_2_i = 1.0
    zeta_2_j = 1.24 # Slater exponent for a Helium 1s orbital
    result_2 = evaluate_1s_nuclear_attraction(zeta_2_i, zeta_2_j)
    
    print(f"\n2. Off-Diagonal Case: Let zeta_i = {zeta_2_i}, zeta_j = {zeta_2_j}")
    print(f"   I = 4 * ({zeta_2_i} * {zeta_2_j})**1.5 / ({zeta_2_i} + {zeta_2_j})**2")
    print(f"   I = 4 * ({zeta_2_i * zeta_2_j})**1.5 / ({zeta_2_i + zeta_2_j})**2")
    print(f"   I = {4 * (zeta_2_i * zeta_2_j)**1.5:.5f} / { (zeta_2_i + zeta_2_j)**2 } = {result_2:.5f}")

if __name__ == '__main__':
    explain_and_calculate()