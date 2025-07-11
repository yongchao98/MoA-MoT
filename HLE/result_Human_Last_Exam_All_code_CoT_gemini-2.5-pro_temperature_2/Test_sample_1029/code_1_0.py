import sympy
from sympy import pprint

def solve_poynting_vector():
    """
    Calculates the Poynting vector for a moving cylindrical rod
    in an external electric field.
    """
    # Define the symbolic variables
    # r: radial distance from the center
    # R: radius of the rod
    # rho: volume charge density
    # E_ext: magnitude of the external electric field
    # v: speed of the rod
    # mu0: permeability of free space
    # epsilon0: permittivity of free space
    r, R, rho, E_ext, v, mu0, epsilon0 = sympy.symbols('r R rho E v mu_0 epsilon_0', real=True, positive=True)

    # --- Step 1: Poynting Vector Inside the Rod (r < R) ---

    # Electric Field inside (E_r_in * r_hat + E_ext * z_hat)
    # From Gauss's Law: E_r * 2*pi*r*L = (rho * pi*r^2*L) / epsilon0
    E_r_in = (rho * r) / (2 * epsilon0)
    E_z_in = E_ext

    # Magnetic Field inside (B_phi_in * phi_hat)
    # From Ampere's Law: B_phi * 2*pi*r = mu0 * I_enc = mu0 * (rho*v) * (pi*r^2)
    B_phi_in = (mu0 * rho * v * r) / 2

    # Poynting Vector S = (1/mu0) * (E x B)
    # S = (1/mu0) * [ (E_r*r_hat + E_z*z_hat) x (B_phi*phi_hat) ]
    # S = (1/mu0) * [ E_r*B_phi*(r_hat x phi_hat) + E_z*B_phi*(z_hat x phi_hat) ]
    # S = (1/mu0) * [ E_r*B_phi*z_hat - E_z*B_phi*r_hat ]
    
    # Calculate components of S_inside
    S_r_in = - (1 / mu0) * E_z_in * B_phi_in
    S_z_in = (1 / mu0) * E_r_in * B_phi_in

    # --- Step 2: Poynting Vector Outside the Rod (r > R) ---

    # Electric Field outside (E_r_out * r_hat + E_ext * z_hat)
    # From Gauss's Law: E_r * 2*pi*r*L = (rho * pi*R^2*L) / epsilon0
    E_r_out = (rho * R**2) / (2 * epsilon0 * r)
    E_z_out = E_ext
    
    # Magnetic Field outside (B_phi_out * phi_hat)
    # From Ampere's Law: B_phi * 2*pi*r = mu0 * I_total = mu0 * (rho*v) * (pi*R^2)
    B_phi_out = (mu0 * rho * v * R**2) / (2 * r)
    
    # Calculate components of S_outside using the same cross product formula
    S_r_out = - (1 / mu0) * E_z_out * B_phi_out
    S_z_out = (1 / mu0) * E_r_out * B_phi_out

    # --- Step 3: Print the results ---
    
    print("The Poynting vector S is given by S = S_r * r_hat + S_z * z_hat")
    print("-" * 60)
    
    print("Inside the rod (r < R):\n")
    print("S_r = ")
    pprint(S_r_in.simplify())
    print("\nS_z = ")
    pprint(S_z_in.simplify())
    print("\nSo, inside the rod, the Poynting vector S_in is:")
    final_S_in = f"({S_r_in.simplify()}) r_hat + ({S_z_in.simplify()}) z_hat"
    print(final_S_in)
    
    print("-" * 60)
    
    print("Outside the rod (r > R):\n")
    print("S_r = ")
    pprint(S_r_out.simplify())
    print("\nS_z = ")
    pprint(S_z_out.simplify())
    print("\nSo, outside the rod, the Poynting vector S_out is:")
    final_S_out = f"({S_r_out.simplify()}) r_hat + ({S_z_out.simplify()}) z_hat"
    print(final_S_out)

solve_poynting_vector()