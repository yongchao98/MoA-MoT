import math
import scipy.constants

def calculate_joule_heat():
    """
    Calculates the Joule heat dissipated from a shrinking, charge-leaking sphere.

    This solution is based on the assumption that the dissipated Joule heat equals
    the initial electrostatic potential energy of the sphere. This occurs if the
    charge leakage happens completely before the sphere shrinks, resulting in zero
    mechanical work done by the electric field.
    """
    try:
        # Get user input for radius and potential
        a = float(input("Enter the initial radius 'a' of the sphere (in meters): "))
        V = float(input("Enter the initial potential 'V' of the sphere (in volts): "))

        if a < 0 or V < 0:
            print("Radius and potential cannot be negative.")
            return

        # Constants
        epsilon_0 = scipy.constants.epsilon_0  # Permittivity of free space
        pi = math.pi

        # The formula for the initial electrostatic energy (and thus the heat)
        # U = 2 * pi * epsilon_0 * a * V^2
        
        # Calculate the final result
        joule_heat = 2 * pi * epsilon_0 * a * (V**2)

        # Output the explanation and the step-by-step calculation
        print("\n--- Calculation ---")
        print("The formula for the total dissipated Joule heat (H), under the assumption W_mech=0, is:")
        print("H = 2 * \u03c0 * \u03b5\u2080 * a * V\u00b2")
        print("\nSubstituting the given values:")
        # Outputting each number in the final equation
        print(f"H = 2 * {pi:.5f} * {epsilon_0:.5e} F/m * {a} m * ({V} V)\u00b2")
        print(f"H = {joule_heat:.4e} Joules")
        
        # Wrapping the final numerical result as requested, although the symbolic answer is more general.
        print(f"\nFinal Answer: {joule_heat}")


    except ValueError:
        print("Invalid input. Please enter numerical values for radius and potential.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    calculate_joule_heat()
