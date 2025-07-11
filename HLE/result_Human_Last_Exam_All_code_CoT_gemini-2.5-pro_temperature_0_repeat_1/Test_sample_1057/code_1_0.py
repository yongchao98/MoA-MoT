import math

def calculate_joule_heat():
    """
    Calculates the Joule heat dissipated by a shrinking, leaking charged sphere.

    This calculation is based on the principle that the total Joule heat dissipated
    is equal to the initial electrostatic potential energy stored in the sphere.
    """

    # --- Constants ---
    # Permittivity of free space in Farads per meter (F/m)
    epsilon_0 = 8.854187817e-12
    pi = math.pi

    # --- Variables (using example values) ---
    # Initial radius of the sphere in meters
    a = 0.1
    # Initial potential of the sphere in Volts
    V = 10000

    # --- Calculation ---
    # The initial stored electrostatic energy U is given by:
    # U = (1/2) * C * V^2
    # where the capacitance C of a sphere is C = 4 * pi * epsilon_0 * a.
    # So, U = (1/2) * (4 * pi * epsilon_0 * a) * V^2 = 2 * pi * epsilon_0 * a * V^2

    joule_heat = 2 * pi * epsilon_0 * a * V**2

    # --- Output ---
    print("The total Joule heat dissipated is equal to the initial electrostatic energy.")
    print("The formula is: Q = 2 * pi * epsilon_0 * a * V^2")
    print("\nUsing the example values:")
    print(f"  Radius (a) = {a} m")
    print(f"  Potential (V) = {V} V")
    print("\nThe final equation with numbers is:")
    # The prompt requires printing each number in the final equation.
    print(f"Q = {2} * {pi} * {epsilon_0} * {a} * {V**2}")

    print(f"\nCalculated Joule Heat (Q) = {joule_heat:.6f} Joules")

if __name__ == "__main__":
    calculate_joule_heat()
