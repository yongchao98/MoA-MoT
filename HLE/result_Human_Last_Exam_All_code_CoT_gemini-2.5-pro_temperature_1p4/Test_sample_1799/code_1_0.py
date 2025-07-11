def solve_relativity_postulate_question():
    """
    Demonstrates why the 2nd postulate of special relativity is not superfluous
    by comparing Galilean and Lorentz velocity addition.
    """
    # The speed of light in a vacuum (meters per second) is a universal constant.
    c = 299792458

    # --- Scenario Setup ---
    # Imagine a spaceship moving away from a stationary observer (e.g., on Earth)
    # at a very high speed: 50% of the speed of light.
    v_spaceship = 0.5 * c

    # The spaceship fires a laser beam in its forward direction. From the perspective
    # of anyone on the spaceship, this light beam travels at speed c.
    v_light_in_ship_frame = c

    print("--- 1. Classical Physics (Galilean Transformation) ---")
    print("This world obeys Postulate 1, but NOT Postulate 2.")
    print("An observer on Earth measures the light's speed by simple addition.")
    print("\nEquation: v_observed = v_spaceship + v_light_in_ship_frame")

    # Perform the Galilean calculation
    v_observed_galilean = v_spaceship + v_light_in_ship_frame

    print(f"Calculation: v_observed = {v_spaceship:.0f} + {v_light_in_ship_frame:.0f}")
    print(f"Result: The observed speed of light would be {v_observed_galilean:.0f} m/s.")
    print(f"This is {v_observed_galilean / c:.1f} times the speed of light, which violates Postulate 2.\n")
    print("-" * 50)

    print("\n--- 2. Special Relativity (Lorentz Transformation) ---")
    print("This world obeys BOTH Postulate 1 and Postulate 2.")
    print("The observer on Earth must use the relativistic formula for velocity addition.")
    print("\nEquation: v_observed = (v_spaceship + v_light_in_ship_frame) / (1 + (v_spaceship * v_light_in_ship_frame) / c^2)")

    # Perform the Lorentz calculation
    numerator = v_spaceship + v_light_in_ship_frame
    denominator = 1 + (v_spaceship * v_light_in_ship_frame) / (c**2)
    v_observed_lorentz = numerator / denominator
    
    # Print the equation with all the numbers
    print(f"Calculation: v_observed = ({v_spaceship:.0f} + {v_light_in_ship_frame:.0f}) / (1 + ({v_spaceship:.0f} * {v_light_in_ship_frame:.0f}) / {c**2:.0f})")
    print(f"Result: The observed speed of light is {v_observed_lorentz:.0f} m/s.")
    print(f"This is exactly equal to c, which is consistent with Postulate 2.")


solve_relativity_postulate_question()