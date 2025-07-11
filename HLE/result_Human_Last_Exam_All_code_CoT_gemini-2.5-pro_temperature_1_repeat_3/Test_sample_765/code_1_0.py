import math
import random

def solve_area():
    """
    Calculates the area of the specified region R using a Monte Carlo method.
    """
    # Plan steps as described in the text.
    # We will use a numerical simulation to estimate the area.

    # Parameters
    num_samples = 10_000_000
    square_side = 6
    domain_area = square_side**2

    # Step 1: Find all integer pairs (m, n) in the domain [0,5]x[0,5]
    # for which sqrt(m^2 + n^2) is an integer C.
    valid_squares = {}
    for m in range(int(square_side)):
        for n in range(int(square_side)):
            c = math.sqrt(m**2 + n**2)
            # Use math.isclose for robust floating-point comparison
            if math.isclose(c, round(c)):
                valid_squares[(m, n)] = round(c)

    # Step 2: Run the Monte Carlo simulation
    points_in_R = 0
    for _ in range(num_samples):
        # Generate a random point in the 6x6 square
        a = random.uniform(0, square_side)
        b = random.uniform(0, square_side)

        # Determine the unit square the point belongs to
        m, n = math.floor(a), math.floor(b)

        # Check if this unit square is a valid one
        if (m, n) in valid_squares:
            C = valid_squares[(m, n)]
            # Check the simplified condition: a^2 + b^2 < (C+1)^2
            if a**2 + b**2 < (C + 1)**2:
                points_in_R += 1

    # Step 3: Calculate the final area
    estimated_area = (points_in_R / num_samples) * domain_area

    # Output the results as requested
    print("Method: Monte Carlo Simulation")
    print(f"1. Total points sampled in the {square_side}x{square_side} square: {num_samples}")
    print(f"2. Points found to be in region R: {points_in_R}")
    print(f"3. Total area of the square: {domain_area}")
    print("\nFinal equation for the area:")
    # The final equation with all numbers
    print(f"Area = ({points_in_R} / {num_samples}) * {domain_area} = {estimated_area:.2f}")

solve_area()