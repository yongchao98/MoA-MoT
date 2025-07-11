def compute_poynting_vector():
    """
    This script symbolically computes the Poynting vector S for a moving,
    charged cylindrical rod in a uniform external electric field.
    The Poynting vector is given by S = (1/mu_0) * (E x B).
    
    The script follows these steps:
    1. Define the components of the total Electric Field (E) and the Magnetic Field (B).
    2. Calculate the components of the cross product (E x B).
    3. Calculate the components of the Poynting vector S by dividing by mu_0.
    4. Print the final result in a clear format.
    """
    
    print("Computing the Poynting vector S = (1/mu_0) * (E x B)")
    print("All calculations are for a point inside the rod (r <= R).")
    print("-" * 50)

    # Step 1: Define the Electric and Magnetic Fields in cylindrical coordinates (r, phi, z).
    # The components are represented as strings for symbolic manipulation.
    print("Step 1: Defining the Electric and Magnetic Fields")
    
    # Total E-field E = E_external + E_rod
    E_vec = {
        "r": "rho * r / (2 * epsilon_0)",
        "phi": "0",
        "z": "E"
    }
    print("Total Electric Field E has components:")
    print(f"  E_r   = {E_vec['r']}")
    print(f"  E_phi = {E_vec['phi']}")
    print(f"  E_z   = {E_vec['z']}\n")

    # B-field from moving charge (J = rho * v * z_hat) via Ampere's Law.
    B_vec = {
        "r": "0",
        "phi": "mu_0 * rho * v * r / 2",
        "z": "0"
    }
    print("Magnetic Field B has components:")
    print(f"  B_r   = {B_vec['r']}")
    print(f"  B_phi = {B_vec['phi']}")
    print(f"  B_z   = {B_vec['z']}")
    print("-" * 50)

    # Step 2: Compute the cross product (E x B).
    # The cross product in cylindrical coordinates is:
    # (A x B)_r   = A_phi*B_z - A_z*B_phi
    # (A x B)_phi = A_z*B_r - A_r*B_z
    # (A x B)_z   = A_r*B_phi - A_phi*B_r
    print("Step 2: Computing the cross product (E x B)")
    
    # After simplification (multiplying by zero terms are removed):
    ExB_r = f"-({E_vec['z']}) * ({B_vec['phi']}) = -E * mu_0 * rho * v * r / 2"
    ExB_phi = "0"
    ExB_z = f"({E_vec['r']}) * ({B_vec['phi']}) = mu_0 * rho**2 * v * r**2 / (4 * epsilon_0)"
    
    print("The components of (E x B) are:")
    print(f"  (E x B)_r   = {ExB_r}")
    print(f"  (E x B)_phi = {ExB_phi}")
    print(f"  (E x B)_z   = {ExB_z}")
    print("-" * 50)

    # Step 3: Compute the Poynting Vector S = (1/mu_0) * (E x B).
    print("Step 3: Computing the final Poynting Vector S")

    # Divide the components of (E x B) by mu_0.
    S_r = "-E * rho * v * r / 2"
    S_phi = "0"
    S_z = "rho**2 * v * r**2 / (4 * epsilon_0)"

    print("The final components of the Poynting vector S are:")
    print(f"  S_r   = {S_r}")
    print(f"  S_phi = {S_phi}")
    print(f"  S_z   = {S_z}\n")

    # Final Output
    print("=" * 50)
    print("FINAL POYNTING VECTOR EQUATION")
    print("=" * 50)
    print("S = (S_r) r_hat + (S_phi) phi_hat + (S_z) z_hat\n")
    print("Substituting the components and showing each number:")
    print(f"S(r) = [ - (1/2) * E * rho * v * r ] r_hat + [ (1/4) * (rho**2 * v * r**2) / epsilon_0 ] z_hat")
    print("=" * 50)

if __name__ == '__main__':
    compute_poynting_vector()