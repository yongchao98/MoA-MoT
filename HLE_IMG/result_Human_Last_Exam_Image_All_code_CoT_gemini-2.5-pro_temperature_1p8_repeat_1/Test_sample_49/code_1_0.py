import sympy as sp

def solve_cutoff_frequency():
    """
    This function symbolically derives and prints the cutoff frequency for the given ladder network.
    """
    # Define symbols for the components and variables
    r, C = sp.symbols('r C', real=True, positive=True)
    Vt = sp.Symbol('Vt') # Test voltage
    k = sp.Symbol('k', integer=True) # Node index

    print("Step 1: Define the cutoff frequency formula")
    print("The cutoff angular frequency, omega_c, is given by:")
    print("omega_c = 1 / (R_eq * C)")
    print("where R_eq is the Thevenin equivalent resistance seen by the capacitor C.")
    print("-" * 40)

    print("Step 2: Find the Thevenin Equivalent Resistance R_eq")
    print("We analyze the infinite resistive ladder with the input (c0) grounded and a test voltage Vt applied at node a0.")
    print("The node voltages V_ck and V_dk can be found by solving the network's difference equations.")
    print("\nThe general solution for the node voltages that decay to zero at infinity is:")
    print("V_ck = A * (2 - sqrt(3))^k + B")
    print("V_dk = -A * (2 - sqrt(3))^k + B")
    print("\nApplying boundary conditions V_c0 = 0 and V_d0 = Vt (at k=0) yields:")
    # A + B = 0  => B = -A
    # -A + B = Vt => -2A = Vt => A = -Vt/2
    # B = Vt/2
    A_val = -Vt / 2
    B_val = Vt / 2
    print(f"A = {A_val}, B = {B_val}")

    # Voltage expressions
    lambda_decay = 2 - sp.sqrt(3)
    V_dk_expr = -A_val * (lambda_decay**k) + B_val

    # Calculate V_d0 and V_d1 to find the current It
    V_d0 = V_dk_expr.subs(k, 0)
    V_d1 = V_dk_expr.subs(k, 1)

    print("\nThe voltage at node d1 is:")
    print(f"V_d1 = {sp.simplify(V_d1)}")
    
    print("\nThe current It from the test source is It = (V_d0 - V_d1) / r:")
    It_expr = (V_d0 - V_d1) / r
    print(f"It = {sp.simplify(It_expr)}")

    print("\nThe equivalent resistance R_eq = Vt / It is:")
    R_eq_expr = sp.simplify(Vt / It_expr)
    print(f"R_eq = {R_eq_expr}")
    print("-" * 40)
    
    print("Step 3: Determine the Cutoff Frequency")
    print("Substituting R_eq into the formula omega_c = 1 / (R_eq * C):")

    # Deconstruct the expression for clean printing
    coeff_term = R_eq_expr / r
    parts = sp.Add.make_args(coeff_term)
    num_1 = parts[0]
    num_sqrt = parts[1]
    radicand = num_sqrt.args[0]
    
    print("\nThe final equation for the cutoff angular frequency is:")
    final_equation = f"omega_c = 1 / (({num_1} + sqrt({radicand})) * r * C)"
    print(final_equation)

    print("\nThe numbers in the final equation's coefficient are:")
    print(f"1) The integer term: {num_1}")
    print(f"2) The number under the square root: {radicand}")

solve_cutoff_frequency()