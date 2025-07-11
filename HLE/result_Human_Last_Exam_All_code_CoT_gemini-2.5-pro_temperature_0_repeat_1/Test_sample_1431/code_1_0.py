def solve_rope_challenge():
    """
    Calculates the force required to lift a rope under specific conditions.

    The problem asks for the force F at the moment a rope of length l and mass m,
    being lifted from one end, has its other end just leaving the ground. At this
    moment, the speed of the end being lifted is v.

    The formula derived from the work-energy principle, including dissipated energy, is:
    F = m * g + (2/3) * m * v^2 / l
    """
    try:
        # Get inputs from the user
        m_str = input("Enter the mass of the rope (m) in kg: ")
        l_str = input("Enter the length of the rope (l) in meters: ")
        v_str = input("Enter the final speed of the top end (v) in m/s: ")

        # Convert inputs to float
        m = float(m_str)
        l = float(l_str)
        v = float(v_str)

        if m <= 0 or l <= 0 or v < 0:
            print("Mass, length must be positive, and speed must be non-negative.")
            return

        # Gravitational constant
        g = 9.8  # m/s^2

        # Calculate the two components of the force
        force_gravity = m * g
        force_kinetic = (2/3) * m * v**2 / l

        # Calculate the total force
        total_force = force_gravity + force_kinetic

        # Print the explanation and the final equation with values
        print("\nThe formula for the force F is: F = m * g + (2/3) * m * v^2 / l")
        print("\nSubstituting the given values:")
        print(f"F = {m} * {g} + (2/3) * {m} * {v}^2 / {l}")
        print(f"F = {force_gravity:.2f} + {force_kinetic:.2f}")
        print(f"The total mystical force F required is: {total_force:.2f} Newtons")

    except ValueError:
        print("Invalid input. Please enter numerical values.")
    except ZeroDivisionError:
        print("Error: Length (l) cannot be zero.")

# Execute the function to solve the challenge
solve_rope_challenge()