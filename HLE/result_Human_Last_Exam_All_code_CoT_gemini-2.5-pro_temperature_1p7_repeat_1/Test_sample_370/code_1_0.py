import sympy

def calculate_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering
    in the high-energy limit for the given Lagrangian, and prints the derivation.
    """

    # Define symbols for the calculation
    g, E, pi = sympy.symbols('g E pi')
    
    print("Theoretical Plan:")
    print("We will calculate the total cross-section for the scattering of two identical fermions: psi + psi -> psi + psi.")
    print("The interaction is given by the Lagrangian term: L_int = -g * bar(psi) * gamma_5 * psi * phi.")
    print("The calculation will be performed in the high-energy limit (E >> m, M).")
    print("The steps are as follows:")
    print("1. Calculate the spin-summed squared matrix elements for the t-channel, u-channel, and their interference term.")
    print("2. Combine these to find the total spin-summed squared matrix element, Sum[|M|^2].")
    print("3. Use the formula for the differential cross-section for identical particles in the center-of-mass (CM) frame.")
    print("4. Integrate over the solid angle to obtain the total cross section, sigma.\n")

    print("Step-by-step Calculation:\n")

    # --- Step 1 & 2: Calculate Sum[|M|^2] ---
    print("Steps 1 & 2: Calculating the total spin-summed squared matrix element Sum[|M|^2]")
    
    # In the high-energy limit, the individual spin-summed squared amplitudes simplify significantly.
    sum_M_t_sq = 4 * g**4
    sum_M_u_sq = 4 * g**4
    
    print(f"The t-channel contribution is: Sum[|M_t|^2] = {sum_M_t_sq}")
    print(f"The u-channel contribution is: Sum[|M_u|^2] = {sum_M_u_sq}")
    
    # The interference term also simplifies in this limit.
    # Note: The asterisk (*) denotes complex conjugation.
    sum_Mt_Mu_star = -2 * g**4
    print(f"The interference term is: Sum[M_t * M_u*] = {sum_Mt_Mu_star}")

    # For identical fermions, the total amplitude is antisymmetrized: M = M_t - M_u.
    total_M_sq_sum = sum_M_t_sq + sum_M_u_sq - 2 * sum_Mt_Mu_star
    print("The total amplitude must be antisymmetrized: M = M_t - M_u.")
    print(f"Sum[|M|^2] = Sum[|M_t|^2] + Sum[|M_u|^2] - 2*Re(Sum[M_t*M_u*])")
    print(f"            = {sum_M_t_sq} + {sum_M_u_sq} - 2*({sum_Mt_Mu_star})")
    print(f"            = {total_M_sq_sum}\n")

    # --- Step 3: Differential Cross Section ---
    print("Step 3: Calculating the differential cross-section d_sigma / d_Omega")
    # In the CM frame, the Mandelstam variable s is (p1+p2)^2 = (2E)^2 = 4E^2.
    s = 4 * E**2
    # The formula for identical final state particles includes a statistical factor of 1/2.
    # d_sigma/d_Omega = (1/2) * (1 / (64*pi^2*s)) * Sum[|M|^2]
    dsigma_dOmega = sympy.Rational(1, 2) * (1 / (64 * pi**2 * s)) * total_M_sq_sum
    print(f"The formula is: d_sigma/d_Omega = (1/2) * (1 / (64*pi^2*s)) * Sum[|M|^2]")
    print(f"With s = 4*E^2 and Sum[|M|^2] = {total_M_sq_sum}, this becomes:")
    print(f"d_sigma/d_Omega = {sympy.simplify(dsigma_dOmega)}\n")

    # --- Step 4: Total Cross Section ---
    print("Step 4: Integrating to find the total cross-section sigma")
    # The differential cross-section is isotropic (independent of angle), so we integrate over the full solid angle (4*pi).
    sigma = 4 * pi * dsigma_dOmega
    sigma_simplified = sympy.simplify(sigma)
    print("sigma = Integral(d_sigma/d_Omega d_Omega) = 4*pi * (d_sigma/d_Omega)")
    print(f"sigma = 4*pi * ({sympy.simplify(dsigma_dOmega)})")
    
    # --- Final Answer ---
    print("\n-----------------------------------------------------")
    print("Final Result:")
    print("The total cross section for the psi-psi scattering is:")
    print(f"sigma = {sigma_simplified}")
    print("\nBreaking down the numbers in the final equation:")
    
    print(f"The numerator has a numerical factor of 3.")
    print(f"The coupling constant 'g' in the numerator is raised to the power of 4.")
    print(f"The denominator has a numerical factor of 32.")
    print(f"The denominator contains 'pi' raised to the power of 1.")
    print(f"The denominator contains the energy 'E' raised to the power of 2.")

if __name__ == '__main__':
    # To run this code, you need to have the sympy library installed (`pip install sympy`).
    calculate_cross_section()
