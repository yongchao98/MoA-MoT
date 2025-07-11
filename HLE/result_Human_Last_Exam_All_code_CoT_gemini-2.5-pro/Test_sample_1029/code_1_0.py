import sympy
from sympy import symbols, pprint, Rational

def compute_poynting_vector():
    """
    This function symbolically computes the Poynting vector at the surface of
    a moving, charged cylindrical rod in an external electric field.
    """
    # Step 1: Define symbolic variables for the physical quantities.
    # R: radius of the rod
    # rho: volume charge density
    # E: external electric field magnitude
    # v: velocity of the rod
    # mu0: permeability of free space
    # eps0: permittivity of free space
    R, rho, E, v, mu0, eps0 = symbols('R ρ E v μ₀ ε₀', real=True, positive=True)

    print("--- Derivation Steps ---")
    
    # Step 2: Define the total Electric Field (E_vec) at the surface (r=R).
    # It has a radial component (from the rod's charge) and an axial component (external field).
    # E_vec = E_r * r̂ + E_z * ẑ
    E_r = (rho * R) / (2 * eps0)
    E_z = E
    print("\n1. Total Electric Field E_vec = E_r * r̂ + E_z * ẑ")
    print("   Radial component (from Gauss's Law), E_r:")
    pprint(E_r)
    print("   Axial component (external field), E_z:")
    pprint(E_z)

    # Step 3: Define the Magnetic Field (B_vec) at the surface (r=R).
    # The moving charge creates a current, generating an azimuthal magnetic field.
    # B_vec = B_phi * φ̂
    B_phi = (mu0 * rho * v * R) / 2
    print("\n2. Magnetic Field B_vec = B_phi * φ̂")
    print("   Azimuthal component (from Ampere's Law), B_phi:")
    pprint(B_phi)

    # Step 4: Compute the cross product E_vec x B_vec in cylindrical coordinates.
    # For A = A_r*r̂ + A_phi*φ̂ + A_z*ẑ and B = B_r*r̂ + B_phi*φ̂ + B_z*ẑ:
    # (A x B)_r = A_phi * B_z - A_z * B_phi
    # (A x B)_phi = A_z * B_r - A_r * B_z
    # (A x B)_z = A_r * B_phi - A_phi * B_r
    # In our case, E_phi = 0, B_r = 0, B_z = 0.
    ExB_r = 0 - E_z * B_phi
    ExB_phi = 0 - 0
    ExB_z = E_r * B_phi - 0
    
    print("\n3. Cross Product E_vec x B_vec")
    print("   Radial component:")
    pprint(ExB_r)
    print("   Axial component:")
    pprint(ExB_z)

    # Step 5: Compute the Poynting vector S_vec = (1/μ₀) * (E_vec x B_vec).
    S_r = sympy.simplify(ExB_r / mu0)
    S_phi = sympy.simplify(ExB_phi / mu0)
    S_z = sympy.simplify(ExB_z / mu0)

    print("\n--- Final Poynting Vector ---")
    print("The Poynting vector S_vec is computed as S_vec = (1/μ₀) * (E_vec x B_vec).")
    print("It has a radial component S_r and an axial component S_z.\n")
    
    # Print the equation for the radial component and its numbers
    print("The final equation for the radial component is: S_r = - (E * ρ * v * R) / 2")
    print("The numbers in this part of the equation are:")
    coeff_r = S_r.as_coeff_Mul()[0]
    print(f"Coefficient: {int(coeff_r.p)} / {int(coeff_r.q)}") # Prints -1 / 2

    # Print the equation for the axial component and its numbers
    print("\nThe final equation for the axial component is: S_z = (ρ² * v * R²) / (4 * ε₀)")
    print("The numbers in this part of the equation are:")
    coeff_z = S_z.as_coeff_Mul()[0]
    print(f"Coefficient: {int(coeff_z.p)} / {int(coeff_z.q)}") # Prints 1 / 4

# Execute the function
compute_poynting_vector()