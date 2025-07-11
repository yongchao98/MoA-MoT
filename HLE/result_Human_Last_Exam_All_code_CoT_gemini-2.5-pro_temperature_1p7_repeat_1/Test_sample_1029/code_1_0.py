import sympy
from sympy.vector import CoordSys3D

def compute_poynting_vector():
    """
    This function symbolically computes the Poynting vector for a moving,
    charged cylindrical rod in an external electric field.
    """

    # 1. Define coordinate system and symbolic variables
    # Using standard cylindrical coordinates (r, phi, z)
    # In sympy.vector, these are C.r, C.theta, C.z with base vectors C.i, C.j, C.k
    # where C.i is r_hat, C.j is phi_hat, C.k is z_hat
    C = CoordSys3D('C', transformation='cylindrical', variable_names=("r", "phi", "z"))
    r, phi, z = C.r, C.phi, C.z
    r_hat, phi_hat, z_hat = C.i, C.j, C.k

    # Define physical constants and variables
    R, rho, E, v, epsilon_0, mu_0 = sympy.symbols('R, rho, E, v, epsilon_0, mu_0', positive=True)

    # --- Step 2: Calculate the Total Electric Field ---

    # External field
    E_ext = E * z_hat

    # Field from the charged rod (using Gauss's Law)
    # Inside the rod (r < R)
    E_rod_inside = (rho * r / (2 * epsilon_0)) * r_hat
    # Outside the rod (r > R)
    E_rod_outside = (rho * R**2 / (2 * epsilon_0 * r)) * r_hat

    # Total electric field
    E_total_inside = E_rod_inside + E_ext
    E_total_outside = E_rod_outside + E_ext

    # --- Step 3: Calculate the Magnetic Field ---

    # From moving charges (using Ampere's Law)
    # B field is zero outside the current source in this calculation method if J=0
    # Inside the rod (r < R)
    B_inside = (mu_0 * rho * v * r / 2) * phi_hat
    # Outside the rod (r > R)
    B_outside = (mu_0 * rho * v * R**2 / (2 * r)) * phi_hat

    # --- Step 4: Compute the Poynting Vector S = (1/mu_0) * (E x B) ---

    # For inside the rod (r < R)
    S_inside = (1 / mu_0) * E_total_inside.cross(B_inside)
    S_inside_simplified = sympy.simplify(S_inside)

    # For outside the rod (r > R)
    S_outside = (1 / mu_0) * E_total_outside.cross(B_outside)
    S_outside_simplified = sympy.simplify(S_outside)
    
    # Extract components for clear printing
    # Reminder: C.i = r_hat, C.k = z_hat
    Sr_inside = S_inside_simplified.dot(r_hat)
    Sz_inside = S_inside_simplified.dot(z_hat)
    
    Sr_outside = S_outside_simplified.dot(r_hat)
    Sz_outside = S_outside_simplified.dot(z_hat)

    # Print the results
    print("The Poynting vector S has two components: a radial component (S_r) and an axial component (S_z).")
    print("S = S_r * r_hat + S_z * z_hat\n")
    
    print("--- For r < R (Inside the rod) ---")
    print(f"S_r = {Sr_inside}")
    print(f"S_z = {Sz_inside}")
    print("\nIn vector form:")
    print(f"S(r<R) = ({Sr_inside})*r_hat + ({Sz_inside})*z_hat")
    
    print("\n" + "="*40 + "\n")
    
    print("--- For r > R (Outside the rod) ---")
    print(f"S_r = {Sr_outside}")
    print(f"S_z = {Sz_outside}")
    print("\nIn vector form:")
    print(f"S(r>R) = ({Sr_outside})*r_hat + ({Sz_outside})*z_hat")


compute_poynting_vector()