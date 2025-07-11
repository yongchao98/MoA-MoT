import math

def galilean_velocity_addition(u, v):
    """
    Calculates velocity addition using the old Galilean formula.
    u: velocity of object relative to moving frame
    v: velocity of moving frame relative to stationary frame
    """
    return u + v

def relativistic_velocity_addition(u, v, c):
    """
    Calculates velocity addition using the correct relativistic formula
    derived from Einstein's postulates.
    u: velocity of object relative to moving frame
    v: velocity of moving frame relative to stationary frame
    c: the speed of light
    """
    numerator = u + v
    denominator = 1 + (u * v) / (c**2)
    return numerator / denominator

def main():
    # The speed of light in a vacuum (meters per second)
    c = 299792458

    # Imagine a spaceship traveling at 90% of the speed of light
    spaceship_velocity = 0.9 * c

    # The spaceship turns on its headlights. The speed of the light beam
    # it emits is 'c' relative to the spaceship.
    light_velocity_relative_to_ship = c

    print("A spaceship is traveling at a velocity v relative to a stationary observer.")
    print("The spaceship emits a beam of light in its direction of travel.")
    print("What is the speed of that light beam as measured by the stationary observer?\n")

    # --- Calculation using the old (incorrect) Galilean physics ---
    print("--- 1. According to Galilean Relativity (Pre-Einstein) ---")
    galilean_result = galilean_velocity_addition(light_velocity_relative_to_ship, spaceship_velocity)

    print(f"The formula is: u' = u + v")
    print(f"The calculation is: u' = {light_velocity_relative_to_ship} + {spaceship_velocity}")
    print(f"Predicted speed of light: {galilean_result:.0f} m/s")
    print(f"This result is {galilean_result/c:.1f} times the speed of light, which violates the 2nd postulate.\n")

    # --- Calculation using Einstein's Special Relativity ---
    print("--- 2. According to Special Relativity (Einstein) ---")
    relativistic_result = relativistic_velocity_addition(light_velocity_relative_to_ship, spaceship_velocity, c)
    
    print(f"The formula is: u' = (u + v) / (1 + (u*v)/c^2)")
    numerator_str = f"({light_velocity_relative_to_ship} + {spaceship_velocity})"
    denominator_str = f"(1 + ({light_velocity_relative_to_ship} * {spaceship_velocity}) / {c}^2)"
    print(f"The calculation is: u' = {numerator_str} / {denominator_str}")
    print(f"Predicted speed of light: {relativistic_result:.0f} m/s")
    print("This result is exactly the speed of light, which is consistent with the 2nd postulate.")


if __name__ == "__main__":
    main()
