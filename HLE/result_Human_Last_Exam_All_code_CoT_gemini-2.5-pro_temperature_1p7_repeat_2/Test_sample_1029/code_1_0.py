import sympy as sp

def solve_poynting_vector():
    """
    This script symbolically computes the Poynting vector for a moving, charged cylindrical rod
    in an external electric field.
    """
    print("Deriving the Poynting Vector S = (1/mu_0) * (E_total x B)")
    print("-" * 60)

    # --- Step 1: Define symbolic variables ---
    # Using sympy to handle the variables symbolically.
    # We will use a cylindrical coordinate system with unit vectors r_hat, phi_hat, z_hat.
    E, rho, v, r, mu_0, epsilon_0 = sp.symbols('E, rho, v, r, mu_0, epsilon_0', positive=True, real=True)
    r_hat = sp.Symbol('r_hat') # Radial unit vector
    phi_hat = sp.Symbol('phi_hat') # Azimuthal unit vector
    z_hat = sp.Symbol('z_hat') # Axial unit vector

    print("Variables defined:")
    print(f"External Electric Field Strength: E")
    print(f"Volume Charge Density: rho")
    print(f"Speed of the rod: v")
    print(f"Radial distance from the axis: r (where r < R)")
    print(f"Permeability of free space: mu_0")
    print(f"Permittivity of free space: epsilon_0")
    print("\n")


    # --- Step 2: Define the Electric and Magnetic Fields ---

    # The total electric field E_total has two components:
    # 1. The external field: E_ext = E * z_hat
    # 2. The field from the rod's charge (via Gauss's Law): E_rod = (rho*r / (2*epsilon_0)) * r_hat
    E_r = rho * r / (2 * epsilon_0)
    E_z = E
    # E_total = E_r * r_hat + E_z * z_hat

    # The magnetic field B is created by the moving charge (current density J = rho*v).
    # From Ampere's Law, B = (mu_0*rho*v*r / 2) * phi_hat
    B_phi = (mu_0 * rho * v * r) / 2
    # B_total = B_phi * phi_hat

    print("Calculating the total Electric and Magnetic fields inside the rod:")
    print(f"Total Electric Field (E_total) = ({sp.pretty(E_r)})*r_hat + ({sp.pretty(E_z)})*z_hat")
    print(f"Magnetic Field (B) = ({sp.pretty(B_phi)})*phi_hat")
    print("\n")

    # --- Step 3: Compute the Poynting Vector ---
    # S = (1/mu_0) * (E_total x B)
    # The cross product in cylindrical coordinates:
    # (E_r*r_hat + E_z*z_hat) x (B_phi*phi_hat)
    # = E_r*B_phi * (r_hat x phi_hat) + E_z*B_phi * (z_hat x phi_hat)
    # Using the identities: r_hat x phi_hat = z_hat and z_hat x phi_hat = -r_hat
    # = (E_r * B_phi)*z_hat - (E_z * B_phi)*r_hat

    # Calculate the components of S
    S_r_component = sp.simplify(-(1/mu_0) * (E_z * B_phi))
    S_z_component = sp.simplify((1/mu_0) * (E_r * B_phi))
    
    print("Calculating the Poynting vector components...")
    
    # Final equation for the radial component (S_r)
    print("S_r = -(1/mu_0) * E_z * B_phi")
    print(f"S_r = -(1/{mu_0}) * ({E}) * ({B_phi})")
    print(f"The radial component of the Poynting vector is:")
    print(f"S_r = {sp.pretty(S_r_component)}")
    print("")

    # Final equation for the axial component (S_z)
    print("S_z = (1/mu_0) * E_r * B_phi")
    print(f"S_z = (1/{mu_0}) * ({E_r}) * ({B_phi})")
    print(f"The axial (z) component of the Poynting vector is:")
    print(f"S_z = {sp.pretty(S_z_component)}")
    print("\n")
    
    # --- Step 4: Final Result ---
    # The final vector is S = S_r*r_hat + S_z*z_hat
    final_vector_str = f"S = ({sp.pretty(S_r_component)})*r_hat + ({sp.pretty(S_z_component)})*z_hat"
    print("-" * 60)
    print("The final Poynting vector is:")
    print(final_vector_str)

    # For the final answer tag, format it as a single string.
    final_answer_str = f"({S_r_component})*r_hat + ({S_z_component})*z_hat"
    return final_answer_str

if __name__ == '__main__':
    final_answer = solve_poynting_vector()
    # The format requested is <<<answer content>>>. The "answer" is the final mathematical expression.
    # print(f"\n<<<{final_answer}>>>")