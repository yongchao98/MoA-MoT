import sympy

def solve_inductance_change():
    """
    This function calculates and prints the expression for the change in mutual inductance.
    """
    # Define the symbolic variables
    mu_0 = sympy.Symbol('μ_0', real=True, positive=True)
    L = sympy.Symbol('L', real=True, positive=True)
    d = sympy.Symbol('d', real=True, positive=True)
    h = sympy.Symbol('h', real=True, positive=True)
    R1 = sympy.Symbol('R1', real=True, positive=True)
    R2 = sympy.Symbol('R2', real=True, positive=True)
    
    # Explain the result for M2
    print("Step 1: Analysis of Mutual Inductance with Concentrators (M2)")
    print("The ideal cylindrical concentrators with radial permeability μ_r -> ∞ and angular permeability μ_θ -> 0")
    print("act as perfect magnetic shields for an internal source. This means the magnetic field outside the shell")
    print("enclosing circuit 1 is zero. Therefore, the flux through circuit 2 is zero.")
    M2 = 0
    print(f"Mutual Inductance with concentrators, M2 = {M2}\n")

    # Calculate M1
    print("Step 2: Analysis of Mutual Inductance for Bare Circuits (M1)")
    print("The mutual inductance M1 between the two bare circuits of length L is given by:")
    # M1 = (mu_0 * L / (2 * sympy.pi)) * sympy.log(d**2 / (d**2 - h**2))
    # An equivalent form is (mu_0 * L / (2 * sympy.pi)) * -log((d**2-h**2)/d**2)
    # The negative form shows that M1 and Delta M have opposite signs.
    # We use the standard positive definition for M1.
    M1_expr = (mu_0 * L / (2 * sympy.pi)) * sympy.log(d**2 / (d**2 - h**2))
    print("M1 = ", end="")
    sympy.pprint(M1_expr, use_unicode=True)
    print("\n")
    
    # Calculate the change in mutual inductance
    print("Step 3: Calculating the Change in Mutual Inductance (ΔM = M2 - M1)")
    delta_M = M2 - M1_expr
    print("ΔM = 0 - M1")
    print("ΔM = ", end="")
    sympy.pprint(delta_M, use_unicode=True)

    # We can also express this in a different but equivalent form
    # delta_M_alt = (mu_0 * L / (2 * sympy.pi)) * sympy.log((d**2 - h**2)/d**2)
    # print("\nOr equivalently:")
    # print("ΔM = ", end="")
    # sympy.pprint(delta_M_alt, use_unicode=True)
    
    # Format the final expression as requested by the prompt.
    final_expr = delta_M
    
    print("\nFinal symbolic expression:")
    print(f"The change in mutual inductance ΔM is given by the equation:")
    # Print each "number" (symbol) in the equation
    # This is a bit of a creative interpretation of the user prompt
    # to print the components of the expression.
    final_str = sympy.sstr(final_expr, full_prec=False)
    # A simple split by operators might not be perfect for complex expressions,
    # but for this one it's illustrative.
    components = final_str.replace(" ","").replace("*"," ").replace("/"," ").replace("("," ").replace(")"," ").replace("-"," - ").split()
    print(f"ΔM = {final_str}")
    # print("The components of this equation are:")
    # for comp in components:
    #     if comp: # check for empty strings
    #         print(comp)

solve_inductance_change()