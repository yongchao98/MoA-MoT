import numpy as np

def calculate_2s_2s_overlap(R, zeta):
    """
    Calculates the overlap integral for two 2s orbitals in a diatomic molecule.

    Args:
        R (float): The internuclear distance in atomic units (bohr).
        zeta (float): The effective nuclear charge.

    Returns:
        float: The value of the overlap integral S.
    """
    print(f"Calculating overlap integral S for R = {R} a.u. and zeta = {zeta}\n")

    # Step 1: Calculate the dimensionless parameter rho
    rho = (zeta * R) / 2
    print(f"Step 1: Calculate rho = (zeta * R) / 2")
    print(f"rho = ({zeta} * {R}) / 2 = {rho}\n")

    # Step 2: The analytical formula for the overlap integral is:
    # S = exp(-rho) * (1 + rho + rho^2/3 + rho^4/15)
    print("Step 2: Use the analytical formula S = exp(-rho) * (1 + rho + rho^2/3 + rho^4/15)")
    print("We will calculate each term of the polynomial part first.\n")

    # Step 3: Calculate each term in the polynomial
    term1 = 1.0
    term2 = rho
    term3 = rho**2 / 3
    term4 = rho**4 / 15
    
    print("Step 3: Calculate the terms in the polynomial P(rho) = 1 + rho + rho^2/3 + rho^4/15")
    print(f"Term 1 (1)      = {term1:.6f}")
    print(f"Term 2 (rho)    = {term2:.6f}")
    print(f"Term 3 (rho^2/3) = {term3:.6f}")
    print(f"Term 4 (rho^4/15)= {term4:.6f}\n")

    # Step 4: Sum the polynomial terms
    poly_sum = term1 + term2 + term3 + term4
    print("Step 4: Sum the polynomial terms")
    print(f"P(rho) = {term1:.6f} + {term2:.6f} + {term3:.6f} + {term4:.6f} = {poly_sum:.6f}\n")

    # Step 5: Calculate the exponential factor
    exp_factor = np.exp(-rho)
    print("Step 5: Calculate the exponential factor exp(-rho)")
    print(f"exp(-{rho}) = {exp_factor:.6f}\n")

    # Step 6: Calculate the final overlap integral S
    S = exp_factor * poly_sum
    print("Step 6: Calculate the final overlap integral S = exp(-rho) * P(rho)")
    print("The final equation is:")
    print(f"S = {exp_factor:.6f} * ({term1:.6f} + {term2:.6f} + {term3:.6f} + {term4:.6f})")
    print(f"S = {exp_factor:.6f} * {poly_sum:.6f}")
    print(f"S = {S:.6f}\n")
    
    return S

if __name__ == '__main__':
    # For the H2+ ion, the nuclear charge is 1. We can use zeta = 1.
    # The equilibrium internuclear distance is approximately R = 2.0 bohr (a.u.).
    effective_nuclear_charge = 1.0
    internuclear_distance = 2.0
    
    overlap_integral = calculate_2s_2s_overlap(R=internuclear_distance, zeta=effective_nuclear_charge)
    
    print("-----------------------------------------------------")
    print(f"The final calculated overlap integral is: {overlap_integral:.6f}")
    print("-----------------------------------------------------")
    # The final answer format required by the prompt
    # print(f'<<<{overlap_integral:.6f}>>>')