import sympy

def solve_joule_heat():
    """
    This function calculates and prints the formula for the Joule heat dissipated
    by a shrinking, charge-leaking sphere.
    """
    # Define the symbols for the variables in the problem
    # a: initial radius of the sphere
    # V: initial potential of the sphere
    # epsilon_0: permittivity of free space
    # pi: the mathematical constant pi
    a, V, epsilon_0, pi = sympy.symbols('a V epsilon_0 pi')

    # Step 1: The initial electrostatic energy of the sphere is U = (1/2) * C * V^2.
    # Step 2: The capacitance C of a sphere of radius 'a' is 4 * pi * epsilon_0 * a.
    C = 4 * pi * epsilon_0 * a

    # Step 3: The initial energy U_initial is the total energy that will be dissipated as Joule heat.
    # The problem states the sphere's radius also decreases, which involves mechanical work (W_mech).
    # The total Joule heat (Q_joule) is Q_joule = U_initial + W_mech.
    # However, the work done depends on the process path (the relative rates of leaking and shrinking).
    # For the problem to have a single, well-defined answer, we must assume a process
    # where the mechanical work done on the field is zero. This happens if the charge
    # leaks away before the sphere shrinks.
    # Therefore, the Joule heat is equal to the initial stored electrostatic energy.
    Q_joule = (sympy.S(1)/2) * C * V**2

    # The final formula for the Joule heat dissipated.
    # We explicitly include the number '2' as requested.
    # The sympy expression simplifies (1/2)*(4*...) to 2*... automatically.
    
    # To ensure the output format is exactly as derived "2 * ...", we can build the string.
    # The simplified expression Q_joule is 2*pi*epsilon_0*a*V**2.
    
    # We will print the components of the final equation.
    coefficient = 2
    term1 = "pi"
    term2 = "epsilon_0"
    term3 = "a"
    term4 = "V**2"
    
    print(f"The total Joule heat dissipated is given by the initial electrostatic energy of the sphere.")
    print(f"The formula for the Joule heat (Q_joule) is:")
    print(f"Q_joule = {coefficient} * {term1} * {term2} * {term3} * {term4}")

solve_joule_heat()