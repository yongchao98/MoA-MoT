import sympy

def solve_cross_section():
    """
    This script calculates the total cross-section for fermion-fermion scattering
    via pseudoscalar meson exchange in the high-energy limit.
    """
    
    # Define symbols for the calculation
    g, E, pi = sympy.symbols('g E pi')
    
    print("This calculation determines the total cross section for psi-psi scattering.")
    print("Assumption: We are in the high-energy limit where particle masses are negligible (m -> 0, M -> 0).\n")
    
    # Step 1: Spin-averaged squared amplitude
    print("Step 1: The spin-averaged squared amplitude |M|^2.")
    print("After applying Feynman rules, accounting for t-channel and u-channel diagrams,")
    print("and performing the spin averaging and summation, the result is surprisingly simple:")
    # In the high-energy limit, |M|^2 = 3 * g^4
    M2_val_num = 3
    M2_val = M2_val_num * g**4
    print(f"|M|^2 = {M2_val_num} * g**4\n")
    
    # Step 2: Differential cross-section in the CM frame
    print("Step 2: The differential cross-section d_sigma / d_Omega.")
    print("The formula in the center-of-mass (CM) frame is: |M|^2 / (64 * pi^2 * s)")
    print("In the CM frame, the Mandelstam variable s = (2*E)^2 = 4*E^2.")
    
    # d_sigma/d_Omega = |M|^2 / (64 * pi^2 * 4 * E^2)
    dsigma_dOmega_num = M2_val_num
    dsigma_dOmega_den = 64 * 4
    
    print(f"d_sigma / d_Omega = ({dsigma_dOmega_num} * g**4) / (64 * pi**2 * 4 * E**2)")
    print(f"d_sigma / d_Omega = {dsigma_dOmega_num} * g**4 / ({dsigma_dOmega_den} * pi**2 * E**2)\n")
    
    # Step 3: Total cross-section
    print("Step 3: The total cross-section sigma.")
    print("For identical particles in the final state, we integrate the differential cross-section")
    print("over a solid angle of 2*pi (one hemisphere).")
    print("sigma = (d_sigma / d_Omega) * 2 * pi")
    
    # sigma = (3 * g^4 / (256 * pi^2 * E^2)) * 2 * pi
    sigma_num = dsigma_dOmega_num * 2
    sigma_den = dsigma_dOmega_den
    
    print(f"sigma = ({dsigma_dOmega_num} * g**4 / ({dsigma_dOmega_den} * pi**2 * E**2)) * 2 * pi")
    print(f"sigma = ({sigma_num} * g**4 * pi) / ({sigma_den} * pi**2 * E**2)")
    
    # Simplify the final expression by cancelling common factors
    final_num = sympy.Integer(sigma_num)
    final_den = sympy.Integer(sigma_den)
    common_divisor = sympy.gcd(final_num, final_den)
    final_num = final_num / common_divisor
    final_den = final_den / common_divisor
    
    print("\nAfter simplifying the numerical factors and pi terms, the final result is:")
    print(f"sigma = {final_num} * g**4 / ({final_den} * pi * E**2)")

solve_cross_section()