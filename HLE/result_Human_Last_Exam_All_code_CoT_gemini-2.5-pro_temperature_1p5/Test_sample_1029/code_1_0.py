import sympy

def compute_poynting_vector():
    """
    This script symbolically computes the Poynting vector S for a moving,
    charged cylindrical rod in an external electric field.
    It calculates the vector components for both inside (r < R) and
    outside (r > R) the rod.
    """
    # Use unicode for prettier output, but fall back to ASCII if not available
    try:
        sympy.init_printing(use_unicode=True)
    except ImportError:
        sympy.init_printing(use_unicode=False)


    # --- Define Symbolic Variables ---
    # mu0: Permeability of free space
    # eps0: Permittivity of free space
    # rho: Volume charge density
    # v: Speed of the rod
    # E_ax: Magnitude of the external axial electric field
    # R: Radius of the rod
    # r: Radial distance from the center of the rod
    mu0, eps0 = sympy.symbols('mu_0 epsilon_0', real=True, positive=True)
    rho, v, E_ax, R = sympy.symbols('rho v E R', real=True, positive=True)
    r = sympy.symbols('r', real=True, positive=True)

    # --- Header for the Output ---
    print("Computation of the Poynting Vector S = (1/mu_0) * (E x B)")
    print("="*60)
    print("A charged rod (radius R, charge density rho) moves at velocity v along its axis.")
    print("An external E-field exists along the axis.\n")

    # --- Region 1: Inside the rod (r < R) ---
    print("--- Inside the rod (r < R) ---")

    # E-field: Sum of external field (axial) and rod's field (radial)
    # From Gauss's Law, E_rod = (rho * r) / (2 * eps0)
    E_r_in = (rho * r) / (2 * eps0)
    E_z_in = E_ax
    
    # B-field: From moving charge (current)
    # From Ampere's Law, B = (mu0 * rho * v * r) / 2
    B_phi_in = (mu0 * rho * v * r) / 2

    # Poynting Vector S = (1/mu0) * (E x B)
    # The cross product (E_z*z_hat + E_r*r_hat) x (B_phi*phi_hat) gives:
    # S = (1/mu0) * [ -E_z*B_phi*r_hat + E_r*B_phi*z_hat ]
    S_r_in = - (E_z_in * B_phi_in) / mu0
    S_z_in = (E_r_in * B_phi_in) / mu0
    
    # Print the resulting components
    print("\nResulting Poynting vector S = S_r * r_hat + S_z * z_hat")
    print("\nRadial component (S_r):")
    sympy.pprint(S_r_in.simplify())
    print("\nAxial component (S_z):")
    sympy.pprint(S_z_in.simplify())
    
    # Create the full vector expression for the final answer
    s_inside_expr = f"({S_r_in.simplify()}) r̂ + ({S_z_in.simplify()}) ẑ"


    print("\n" + "="*60 + "\n")

    # --- Region 2: Outside the rod (r > R) ---
    print("--- Outside the rod (r > R) ---")

    # E-field: Sum of external field (axial) and rod's field (radial)
    # From Gauss's Law, E_rod = (rho * R^2) / (2 * eps0 * r)
    E_r_out = (rho * R**2) / (2 * eps0 * r)
    E_z_out = E_ax
    
    # B-field: From moving charge (current)
    # From Ampere's Law, B = (mu0 * rho * v * R^2) / (2 * r)
    B_phi_out = (mu0 * rho * v * R**2) / (2 * r)

    # Poynting Vector S = (1/mu0) * (E x B)
    S_r_out = - (E_z_out * B_phi_out) / mu0
    S_z_out = (E_r_out * B_phi_out) / mu0

    # Print the resulting components
    print("\nResulting Poynting vector S = S_r * r_hat + S_z * z_hat")
    print("\nRadial component (S_r):")
    sympy.pprint(S_r_out.simplify())
    print("\nAxial component (S_z):")
    sympy.pprint(S_z_out.simplify())

    s_outside_expr = f"({S_r_out.simplify()}) r̂ + ({S_z_out.simplify()}) ẑ"

    # For the final answer format
    print(f"\n\n<<<Inside (r<R): S = {s_inside_expr}\nOutside (r>R): S = {s_outside_expr}>>>")

if __name__ == '__main__':
    compute_poynting_vector()
