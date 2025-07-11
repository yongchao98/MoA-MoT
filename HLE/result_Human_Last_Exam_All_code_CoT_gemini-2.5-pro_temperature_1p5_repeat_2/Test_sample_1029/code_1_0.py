import sympy
from sympy.vector import CoordSys3D

def solve_poynting_vector():
    """
    This script symbolically computes the Poynting vector S for a moving,
    charged cylindrical rod in an external electric field.

    The calculation is performed for two regions:
    1. Inside the rod (r < R)
    2. Outside the rod (r > R)
    """

    # 1. Define a cylindrical coordinate system and symbols
    # We use a standard 3D Cartesian system and will interpret the
    # unit vectors as cylindrical ones: C.i -> r_hat, C.j -> phi_hat, C.k -> z_hat
    C = CoordSys3D('C')
    r_hat, phi_hat, z_hat = C.i, C.j, C.k

    # Define physical constants and variables as symbols
    mu0, epsilon0, rho, v, E_mag, R, r = sympy.symbols('μ₀ ε₀ ρ v E R r', positive=True, real=True)
    
    print("Derivation of the Poynting Vector S = (1/μ₀) * E_total × B\n")

    # --- Calculation for INSIDE the rod (r < R) ---

    print("--- For inside the rod (r < R) ---")

    # 2. Define the Electric and Magnetic Fields
    # External Electric Field (along the axis, z_hat direction)
    E_ext = E_mag * z_hat
    
    # Electric field from the rod's own charge (radial, found with Gauss's Law)
    E_rod_inside = (rho * r / (2 * epsilon0)) * r_hat

    # Total Electric Field is the vector sum
    E_total_inside = E_ext + E_rod_inside
    
    # Magnetic Field from the moving charge (azimuthal, found with Ampere's Law)
    # The current density is J = ρ*v, and the enclosed current is I_enc = J * π*r²
    B_inside = (mu0 * rho * v * r / 2) * phi_hat

    # 3. Compute the Poynting Vector
    # S = (1/μ₀) * E × B
    S_inside = (1 / mu0) * E_total_inside.cross(B_inside)
    S_inside_simplified = sympy.simplify(S_inside)

    # 4. Print the final result for the inside region
    print(f"Total Electric Field: E_total = {E_rod_inside} + {E_ext}")
    print(f"Magnetic Field: B = {B_inside}")
    print("\nResulting Poynting Vector:")
    
    # To present the result in a standard format
    radial_comp_in = S_inside_simplified.dot(r_hat)
    axial_comp_in = S_inside_simplified.dot(z_hat)
    print(f"S_inside = ({radial_comp_in}) r̂ + ({axial_comp_in}) ẑ")
    print("-" * 40)

    # --- Calculation for OUTSIDE the rod (r > R) ---

    print("\n--- For outside the rod (r > R) ---")

    # 2. Define the Electric and Magnetic Fields
    # Electric field from the rod's own charge (Gauss's Law)
    E_rod_outside = (rho * R**2 / (2 * epsilon0 * r)) * r_hat

    # Total Electric Field is the vector sum
    E_total_outside = E_ext + E_rod_outside

    # Magnetic Field from the moving charge (Ampere's Law)
    # Total current is I_total = ρ*v * π*R²
    B_outside = (mu0 * rho * v * R**2 / (2 * r)) * phi_hat

    # 3. Compute the Poynting Vector
    S_outside = (1 / mu0) * E_total_outside.cross(B_outside)
    S_outside_simplified = sympy.simplify(S_outside)

    # 4. Print the final result for the outside region
    print(f"Total Electric Field: E_total = {E_rod_outside} + {E_ext}")
    print(f"Magnetic Field: B = {B_outside}")
    print("\nResulting Poynting Vector:")
    
    radial_comp_out = S_outside_simplified.dot(r_hat)
    axial_comp_out = S_outside_simplified.dot(z_hat)
    print(f"S_outside = ({radial_comp_out}) r̂ + ({axial_comp_out}) ẑ")

if __name__ == '__main__':
    solve_poynting_vector()