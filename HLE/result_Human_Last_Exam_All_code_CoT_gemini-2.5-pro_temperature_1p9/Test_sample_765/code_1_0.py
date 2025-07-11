import math
import random

def calculate_area_of_R():
    """
    Calculates the area of the region R using the Monte Carlo method.

    The region R in the Complex plane is defined by floor(|z|) = |floor(z)|
    within the square defined by the corners (0+0i, 6+6i).

    Let z = x + yi.
    The condition is floor(sqrt(x^2 + y^2)) == sqrt(floor(x)^2 + floor(y)^2).
    """

    num_samples = 10_000_000
    hits = 0
    square_side = 6.0
    total_area = square_side * square_side

    for _ in range(num_samples):
        # Generate a random point z = x + yi within the 6x6 square
        x = random.uniform(0, square_side)
        y = random.uniform(0, square_side)

        # Let n = floor(x) and m = floor(y)
        n = math.floor(x)
        m = math.floor(y)

        # For the equality to hold, |floor(z)| must be an integer.
        # This means n^2 + m^2 must be a perfect square.
        k_squared = n**2 + m**2
        k = math.sqrt(k_squared)

        # Check if k is an integer (using a small tolerance for floating point errors)
        if abs(k - round(k)) < 1e-9:
            # If k is an integer, proceed to check the full condition
            rhs = round(k)
            
            # Calculate the left-hand side of the equation
            mod_z = math.sqrt(x**2 + y**2)
            lhs = math.floor(mod_z)

            # If the condition holds, the point is in the region R
            if lhs == rhs:
                hits += 1

    # Estimate the area
    estimated_area = (hits / num_samples) * total_area
    
    # Print the equation as requested
    print(f"Area = ({hits} / {num_samples}) * {total_area:.2f}")

    # Print the final result
    print(f"The estimated area of R is: {estimated_area:.2f}")
    
    # Return the value for the final answer block
    return estimated_area

# Run the calculation and store the result
final_area = calculate_area_of_R()
# The required format for the final answer
# print(f"<<<{final_area:.2f}>>>")