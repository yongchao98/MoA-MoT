import sympy
from sympy.vector import CoordSys3D

def compute_poynting_vector():
    """
    Computes the Poynting vector inside and outside a moving, charged cylindrical rod
    in an external electric field.
    """
    # Define the symbolic variables for the problem
    # r: radial distance from the axis
    # R: radius of the rod
    # rho: volume charge density
    # v: speed of the rod
    # E: magnitude of the external electric field
    # mu_0: permeability of free space
    # epsilon_0: permittivity of free space
    r, R, rho, v, E = sympy.symbols('r R rho v E', real=True, positive=True)
    mu_0, epsilon_0 = sympy.symbols('mu_0 epsilon_0', real=True, positive=True)

    # Set up a cylindrical coordinate system
    # C.i, C.j, C.k correspond to the unit vectors r_hat, phi_hat, z_hat
    # C.r, C.phi, C.z correspond to the coordinate variables
    C = CoordSys3D('C', transformation='cylindrical', variable_names=("r", "phi", "z"))

    # --- Case 1: Inside the rod (r < R) ---
    print("--- Inside the rod (r < R) ---")

    # Total Electric Field (E_total = E_external + E_rod)
    E_ext = E * C.k  # External field along z-axis
    E_rod_in = (rho * C.r) / (2 * epsilon_0) * C.i # From Gauss's Law
    E_total_in = E_ext + E_rod_in

    # Magnetic Field (from moving charge)
    B_in = (mu_0 * rho * v * C.r) / 2 * C.j # From Ampere's Law

    # Poynting Vector Calculation: S = (1/mu_0) * (E x B)
    S_in = (1 / mu_0) * (E_total_in.cross(B_in))
    S_in_simplified = sympy.simplify(S_in)

    # Output the components of the Poynting vector
    print("The Poynting vector S is expressed in cylindrical coordinates (r_hat, phi_hat, z_hat).")
    print("S = (S_r) r_hat + (S_phi) phi_hat + (S_z) z_hat")
    print(f"S_r = {S_in_simplified.dot(C.i)}")
    print(f"S_phi = {S_in_simplified.dot(C.j)}")
    print(f"S_z = {S_in_simplified.dot(C.k)}")
    print("\n" + "="*40 + "\n")

    # --- Case 2: Outside the rod (r > R) ---
    print("--- Outside the rod (r > R) ---")

    # Total Electric Field (E_total = E_external + E_rod)
    E_rod_out = (rho * R**2) / (2 * epsilon_0 * C.r) * C.i # From Gauss's Law
    E_total_out = E_ext + E_rod_out

    # Magnetic Field (from moving charge)
    B_out = (mu_0 * rho * v * R**2) / (2 * C.r) * C.j # From Ampere's Law

    # Poynting Vector Calculation: S = (1/mu_0) * (E x B)
    S_out = (1 / mu_0) * (E_total_out.cross(B_out))
    S_out_simplified = sympy.simplify(S_out)

    # Output the components of the Poynting vector
    print("The Poynting vector S is expressed in cylindrical coordinates (r_hat, phi_hat, z_hat).")
    print("S = (S_r) r_hat + (S_phi) phi_hat + (S_z) z_hat")
    print(f"S_r = {S_out_simplified.dot(C.i)}")
    print(f"S_phi = {S_out_simplified.dot(C.j)}")
    print(f"S_z = {S_out_simplified.dot(C.k)}")

    # Construct the final answer string for the <<<>>> block
    s_in_str = f"Inside (r<R): S = ({S_in_simplified.dot(C.i)}) r_hat + ({S_in_simplified.dot(C.k)}) z_hat"
    s_out_str = f"Outside (r>R): S = ({S_out_simplified.dot(C.i)}) r_hat + ({S_out_simplified.dot(C.k)}) z_hat"
    final_answer = f"<<<{s_in_str}; {s_out_str}>>>"
    print(f"\n{final_answer}")


if __name__ == '__main__':
    compute_poynting_vector()