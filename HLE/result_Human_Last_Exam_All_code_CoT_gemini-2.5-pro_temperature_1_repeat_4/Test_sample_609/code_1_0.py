import math

def calculate_area_ratio(n):
    """
    Calculates the area ratio between an n-sided polygon formed by extending
    the alternate sides of a 2n-sided polygon.

    The function prints a step-by-step calculation for a given n.
    """
    if n < 3:
        print("Error: The number of sides 'n' must be 3 or greater.")
        return

    # Angle in degrees for the formula 1 / (1 - tan^2(90/n))
    angle_in_degrees = 90.0 / n

    # Convert angle to radians for Python's math functions
    angle_in_radians = math.radians(angle_in_degrees)

    # Calculate tan(angle)^2
    tan_squared = math.tan(angle_in_radians)**2

    # Calculate the final ratio
    ratio = 1 / (1 - tan_squared)

    # Print the results showing the numbers used in the final equation
    print(f"For a starting 2n-sided polygon where n = {n}:")
    print(f"The calculation is based on the formula: 1 / (1 - tan^2(90/n))")
    print("\n--- Calculation Steps ---")
    print(f"1. Plug in n: 1 / (1 - tan^2(90/{n}))")
    print(f"2. Evaluate angle: 1 / (1 - tan^2({angle_in_degrees:.2f}°))")
    print(f"3. Evaluate tan^2 term: 1 / (1 - {tan_squared:.4f})")
    print(f"4. Final Result: {ratio:.4f}")
    
    print("\n--- Final Equation ---")
    # This line shows the final equation with the numbers plugged in as requested.
    print(f"Area Ratio = 1 / (1 - tan^2(90/{n})) = {ratio}")


# Run the calculation for the example case given in the prompt (n=3)
calculate_area_ratio(3)