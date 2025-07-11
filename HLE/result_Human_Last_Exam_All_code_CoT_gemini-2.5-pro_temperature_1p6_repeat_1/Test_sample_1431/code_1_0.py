import sys

def solve_rope_challenge():
    """
    Calculates the force F required to lift a rope of length l and mass m,
    such that its speed is v when the end just leaves the ground.
    """
    try:
        # User inputs for the rope's properties
        m_str = input("Enter the mass of the rope (m) in kilograms: ")
        m = float(m_str)

        l_str = input("Enter the length of the rope (l) in meters: ")
        l = float(l_str)

        v_str = input("Enter the final speed of the rope (v) in meters per second: ")
        v = float(v_str)

        if l <= 0:
            print("Error: Length (l) must be a positive number.", file=sys.stderr)
            return

    except ValueError:
        print("Error: Invalid input. Please enter valid numbers for mass, length, and speed.", file=sys.stderr)
        return

    # Standard gravity constant
    g = 9.8  # m/s^2

    # Calculate the force F using the work-energy theorem
    # F * l - m*g*l/2 = 0.5 * m * v**2
    # F = (0.5 * m * v**2 + 0.5 * m * g * l) / l
    # F = m*v**2/(2*l) + m*g/2
    force_F = (m * v**2) / (2 * l) + (m * g) / 2

    # Print the final equation with the user-provided numbers
    print("\nBased on the work-energy principle, the force F is calculated using the formula:")
    print("F = (m * v^2) / (2 * l) + (m * g) / 2")
    
    print("\nPlugging in your values:")
    # Using format specifiers to ensure all numbers are displayed as floats
    print(f"F = ({m:.1f} * {v:.1f}**2) / (2 * {l:.1f}) + ({m:.1f} * {g:.1f}) / 2")
    
    # Print the final calculated result
    print(f"\nThe calculated mystical force F required is: {force_F:.2f} Newtons")


if __name__ == "__main__":
    solve_rope_challenge()