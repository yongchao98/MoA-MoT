import sys

def solve_energy_shift():
    """
    This function calculates and prints the ground state energy shift
    for two interacting quantum harmonic oscillators.
    """

    # Define the symbols used in the equations as strings
    e = "e"
    m = "m"
    omega_0 = "w_0"
    R = "R"
    hbar = "hbar"
    epsilon_0 = "epsilon_0"
    pi = "pi"
    x1 = "x_1"
    x2 = "x_2"

    print("Step 1: The Dipole-Dipole Interaction Potential (V_int)")
    print("---------------------------------------------------------")
    print(f"The two quantum harmonic oscillators are modeled as electric dipoles.")
    print(f"For a large separation R, the Coulomb interaction is approximated by the dipole-dipole potential.")
    print(f"Assuming the dipoles (of moment p=e*x) are aligned along the line connecting their centers, the potential is:")
    v_int_str = f"- (2 * {e}^2 * {x1} * {x2}) / (4 * {pi} * {epsilon_0} * {R}^3)"
    v_int_simplified = f"- ({e}^2 * {x1} * {x2}) / (2 * {pi} * {epsilon_0} * {R}^3)"
    print(f"V_int = {v_int_str}")
    print(f"This simplifies to the interaction Hamiltonian we will use:")
    print(f"H' = V_int = {v_int_simplified}\n")

    print("Step 2: First-Order Energy Correction")
    print("--------------------------------------")
    print("The first-order energy shift is the expectation value of H' in the unperturbed ground state |0,0>.")
    print("Delta_E^(1) = <0,0| H' |0,0> = - <0|x_1|0> * <0|x_2|0> / (2*pi*epsilon_0*R^3)")
    print("For a quantum harmonic oscillator, the expectation value of position x in any energy eigenstate is 0.")
    print("Therefore, <0|x|0> = 0, and the first-order correction is zero.")
    print("Delta_E^(1) = 0\n")

    print("Step 3: Second-Order Energy Correction")
    print("---------------------------------------")
    print("We must calculate the second-order correction, which is the leading non-zero term.")
    print("The formula is: Delta_E^(2) = sum_{k!=0} |<k| H' |0,0>|^2 / (E_0 - E_k)\n")

    print("Step 4: The Matrix Element <k| H' |0,0>")
    print("-----------------------------------------")
    print("H' is proportional to x_1 * x_2. The position operator x can be written in terms of creation (a_dag) and annihilation (a) operators:")
    print(f"x = sqrt({hbar}/(2*{m}*{omega_0})) * (a + a_dag)")
    print("The operator x only connects states with a change in quantum number by 1 (Delta_n = +/- 1).")
    print("Therefore, the matrix element <n_1,n_2| x_1*x_2 |0,0> is only non-zero for n_1=1 and n_2=1.")
    print("So the sum reduces to a single term with the state |k> = |1,1>.")
    x_10_str = f"sqrt({hbar} / (2 * {m} * {omega_0}))"
    print(f"The specific matrix element required is <1|x|0> = {x_10_str}")
    
    C_factor = f"-({e}^2 / (2 * {pi} * {epsilon_0} * {R}^3))"
    matrix_element = f"{C_factor} * <1|{x1}|0> * <1|{x2}|0>"
    matrix_element_val = f"{C_factor} * ({hbar} / (2 * {m} * {omega_0}))"
    matrix_element_simplified = f"-({e}^2 * {hbar}) / (4 * {pi} * {epsilon_0} * {R}^3 * {m} * {omega_0})"
    print(f"The full matrix element is M = <1,1| H' |0,0> = {matrix_element}")
    print(f"M = {matrix_element_val}")
    print(f"M = {matrix_element_simplified}\n")

    print("Step 5: The Energy Denominator (E_0 - E_k)")
    print("---------------------------------------------")
    e_n = f"(n + 1/2)*{hbar}*{omega_0}"
    print(f"The energy of state |n1, n2> is E_n1,n2 = {e_n.replace('n', 'n1')} + {e_n.replace('n', 'n2')}")
    print(f"Ground state energy E_0 = E_0,0 = {hbar}*{omega_0}")
    print(f"Excited state energy E_k = E_1,1 = 3*{hbar}*{omega_0}")
    energy_denom = f"-2*{hbar}*{omega_0}"
    print(f"The denominator is E_0 - E_k = {hbar}*{omega_0} - 3*{hbar}*{omega_0} = {energy_denom}\n")

    print("Step 6: The Final Energy Shift")
    print("-------------------------------")
    print("Now we assemble the final result: Delta_E = |M|^2 / (E_0 - E_k)")
    
    m_squared_str = f"({e}^4 * {hbar}^2) / (16 * {pi}^2 * {epsilon_0}^2 * {R}^6 * {m}^2 * {omega_0}^2)"
    print(f"|M|^2 = {m_squared_str}")
    
    final_result = f"- ({e}^4 * {hbar}) / (32 * {pi}^2 * {epsilon_0}^2 * {m}^2 * {omega_0}^3 * {R}^6)"
    print(f"\nDelta_E = ( {m_squared_str} ) / ( {energy_denom} )")
    print("\nThe final ground state energy shift is:")
    print(f"Delta_E = {final_result}")

solve_energy_shift()