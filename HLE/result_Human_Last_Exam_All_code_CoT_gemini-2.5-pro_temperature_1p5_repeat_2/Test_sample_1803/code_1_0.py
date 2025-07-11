import sympy

def calculate_energy_shift():
    """
    Calculates and prints the ground state energy shift for two interacting
    quantum harmonic oscillators using second-order perturbation theory.
    """

    # Define the physical constants and variables as symbolic objects
    e, m, omega_0, R, hbar = sympy.symbols('e m omega_0 R hbar', real=True, positive=True)
    pi = sympy.pi

    print("Step 1: Define the dipole-dipole interaction potential V.")
    print("----------------------------------------------------------")
    # For two dipoles p1=e*x1 and p2=e*x2 separated by R along the z-axis,
    # the potential is V = (1/(4*pi*R^3)) * [p1.p2 - 3(p1.z_hat)(p2.z_hat)].
    # The problem specifies using V = q1*q2 / (4*pi*r).
    # This leads to V = (e**2 / (4*pi*R**3)) * [x1*x2 + y1*y2 - 2*z1*z2]
    # We define C as the constant factor in this potential.
    C = (e**2) / (4 * pi * R**3)
    print("The interaction potential V is approximated by the dipole-dipole interaction:")
    print(f"V ≈ C * [x1*x2 + y1*y2 - 2*z1*z2], where C = {C}\n")

    print("Step 2: Calculate the second-order energy shift.")
    print("--------------------------------------------------")
    # The first-order shift <0|V|0> is zero.
    # The second-order shift is Delta_E = sum_{k!=0} |<k|V|0>|^2 / (E0 - Ek).
    # The operator V couples the ground state |00> to states |1i,1j>, where
    # both oscillators are excited. For our V, i must equal j.
    # The energy of the intermediate state |1,1> is E0 + 2*hbar*omega_0.
    E_denominator = -2 * hbar * omega_0
    print("The first-order energy shift is zero.")
    print(f"The energy denominator for the relevant transitions is E0 - E_intermediate = {E_denominator}\n")


    print("Step 3: Calculate the squared matrix elements |<k|V|0>|^2.")
    print("---------------------------------------------------------")
    # The matrix element <1|x|0> for a QHO is sqrt(hbar/(2*m*omega_0)).
    # So, <1i,1i| V_i |0,0> = (coeff_i * C) * (<1|xi|0>)^2
    # where coeff_x=1, coeff_y=1, coeff_z=-2.
    matrix_element_x0_sq = (hbar / (2 * m * omega_0))
    
    # Squared matrix element for the x-x interaction part
    M_x_sq = (C**2) * matrix_element_x0_sq**2
    
    # Squared matrix element for the y-y interaction part
    M_y_sq = (C**2) * matrix_element_x0_sq**2
    
    # Squared matrix element for the z-z interaction part (coefficient is -2)
    M_z_sq = ((-2 * C)**2) * matrix_element_x0_sq**2

    print(f"The squared matrix element for the x-term: |M_x|^2 = {sympy.simplify(M_x_sq)}")
    print(f"The squared matrix element for the y-term: |M_y|^2 = {sympy.simplify(M_y_sq)}")
    print(f"The squared matrix element for the z-term: |M_z|^2 = {sympy.simplify(M_z_sq)}\n")

    print("Step 4: Sum the contributions and derive the final expression.")
    print("--------------------------------------------------------------")
    # Delta_E = ( |M_x|^2 + |M_y|^2 + |M_z|^2 ) / E_denominator
    total_matrix_element_sq = M_x_sq + M_y_sq + M_z_sq
    print(f"The sum of squared matrix elements is 1^2 + 1^2 + (-2)^2 = 6 times the base value:")
    print(f"Sum |M|^2 = {sympy.simplify(total_matrix_element_sq)}")
    
    Delta_E = total_matrix_element_sq / E_denominator
    final_expression = sympy.simplify(Delta_E)

    print("\nThe total ground state energy shift is the sum of these contributions divided by the energy denominator.")
    print("Final Result:")
    # The problem asks to output each number, so we print a formatted string.
    # The simplified result is -3*hbar*e**4 / (64*pi**2*m**2*omega_0**3*R**6)
    # The numbers are 3, 4, 64 = 4*4^2, etc.
    print(f"ΔE₀ = - (3 * ħ * e⁴) / (64 * π² * m² * ω₀³ * R⁶)")
    
if __name__ == '__main__':
    calculate_energy_shift()