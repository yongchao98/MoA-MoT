import math

def calculate_lifting_force():
    """
    Calculates the force required to lift a rope to a certain speed.

    The calculation is based on the work-energy principle:
    Work_done_by_force = Change_in_Kinetic_Energy + Change_in_Potential_Energy
    F * l = (1/2 * m * v^2) + (m * g * l / 2)
    F = m * g / 2 + m * v^2 / (2 * l)
    """
    try:
        # Get user inputs for mass, length, and final speed
        m = float(input("Enter the mass of the rope (m) in kilograms: "))
        l = float(input("Enter the length of the rope (l) in meters: "))
        v = float(input("Enter the final speed of the rope (v) in meters/second: "))

        if m <= 0 or l <= 0 or v < 0:
            print("\nMass, length must be positive, and speed must be non-negative. Please try again.")
            return

        # Standard acceleration due to gravity
        g = 9.8  # m/s^2

        # Calculate the two components of the force
        force_due_to_pe = (m * g) / 2
        force_due_to_ke = (m * v**2) / (2 * l)

        # Calculate the total force
        total_force = force_due_to_pe + force_due_to_ke

        # Print the detailed breakdown of the calculation
        print("\n--- Calculation Breakdown ---")
        print(f"The formula for the force F is: F = (m * g / 2) + (m * v^2 / (2 * l))")
        print("\nGiven values:")
        print(f"Mass (m) = {m:.1f} kg")
        print(f"Length (l) = {l:.1f} m")
        print(f"Final Speed (v) = {v:.1f} m/s")
        print(f"Gravity (g) = {g:.1f} m/s^2")
        print("\nPlugging the numbers into the formula:")
        print(f"F = ({m:.1f} * {g:.1f} / 2) + ({m:.1f} * {v:.1f}**2 / (2 * {l:.1f}))")
        print(f"F = {force_due_to_pe:.2f} N (from potential energy) + {force_due_to_ke:.2f} N (from kinetic energy)")
        print(f"F = {total_force:.2f} N")
        print("\n--- Final Answer ---")
        print(f"The exact mystical force F required is {total_force:.2f} Newtons.")

    except ValueError:
        print("\nInvalid input. Please enter numerical values.")
    except ZeroDivisionError:
        print("\nError: The length of the rope cannot be zero.")

if __name__ == "__main__":
    calculate_lifting_force()