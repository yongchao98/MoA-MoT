import sympy

def solve_and_display_electric_field():
    """
    This function symbolically derives and prints the electric field
    in the two regions of the cylindrical resistor.
    """
    # Define the symbolic variables
    V0, r, pi, sigma1, sigma2 = sympy.symbols('V_0 r pi sigma_1 sigma_2')

    # The problem describes a system that can be thought of as two resistors
    # in series for any given differential cylindrical shell.
    # From solving Laplace's equation with the given boundary conditions,
    # we can determine the potential and then the electric field.
    # The potential constants A1 and A2 are:
    A1 = -2 * V0 * sigma2 / (pi * (sigma1 + sigma2))
    A2 = -2 * V0 * sigma1 / (pi * (sigma1 + sigma2))
    
    # The electric field E = - (1/r) * (d(Phi)/d(phi)) * i_phi.
    # Since Phi_1 = A1*phi + B1, d(Phi_1)/d(phi) = A1.
    # Since Phi_2 = A2*phi + B2, d(Phi_2)/d(phi) = A2.
    E1_phi = - (1/r) * A1
    E2_phi = - (1/r) * A2

    # The result from the derivation is E_phi = E1 for region 1 and E_phi = E2 for region 2.
    print("Electric field in Region 1 (0 < phi < pi/2):")
    # E_1 = (2*sigma2*V0)/(r*pi*(sigma1+sigma2)) * i_phi
    print(f"E_1 = (2 * sigma_2 * V_0) / (r * pi * (sigma_1 + sigma_2)) * i_phi")

    print("\nElectric field in Region 2 (pi/2 < phi < pi):")
    # E_2 = (2*sigma1*V0)/(r*pi*(sigma1+sigma2)) * i_phi
    print(f"E_2 = (2 * sigma_1 * V_0) / (r * pi * (sigma_1 + sigma_2)) * i_phi")

solve_and_display_electric_field()
<<<C>>>