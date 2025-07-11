import math

def display_area_formula():
    """
    This function derives and prints the formula for the area of the triangle T(t).
    """

    # --- 1. Define Constants from the Problem ---
    # Radius of the circle C
    R = 10.0
    # Speed of the triangle's vertices
    v = 1.0
    # Side length of the hexagon H is equal to the radius R
    s = R

    # --- 2. Explain the Derivation Logic ---
    print("Finding the area of the triangle T(t) as a function of time t.")
    print("-" * 60)
    print("Key Principles:")
    print("1. The rotation of the hexagon does not change the triangle's area.")
    print("2. Due to symmetry, the triangle T(t) is always equilateral.")
    print("3. The area of an equilateral triangle with side L is: Area = (sqrt(3)/4) * L^2.")
    print("\nOur goal is to find the side length squared, L(t)^2.")

    # --- 3. Calculate L(t)^2 ---
    # The side length squared L(t)^2 can be found using coordinate geometry.
    # It has been derived that L(t)^2 is a function of the initial side length squared L(0)^2
    # and the distance moved by the vertices (v*t).
    # L(0) is the side length when vertices are at the midpoints, which is 1.5 * s.
    L0_squared = (1.5 * s)**2
    # The time-dependent part of the formula is 3 * (v*t)^2.
    time_factor = 3 * v**2

    print("\nThe side length squared, L(t)^2, is calculated as:")
    print(f"L(t)^2 = (1.5 * s)^2 + 3 * (v*t)^2")
    print(f"L(t)^2 = (1.5 * {s})^2 + {int(time_factor)} * ({int(v)}*t)^2")
    print(f"L(t)^2 = {L0_squared} + {int(time_factor)} * t^2")
    print("-" * 60)

    # --- 4. Present the Final Area Formula ---
    print("\nSubstituting L(t)^2 into the area formula gives the final result.")
    print("The final equation for the area as a function of time is:")
    
    # Using strings to represent the symbolic parts of the equation
    c1_num = "sqrt(3)"
    c1_den = 4
    c2 = L0_squared
    c3 = int(time_factor)

    print(f"\nArea(t) = ({c1_num} / {c1_den}) * ({c2} + {c3} * t^2)\n")


if __name__ == "__main__":
    display_area_formula()