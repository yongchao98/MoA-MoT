import math

def solve_parallelogram_diagonals():
    """
    This function solves the problem of finding the product a*b based on the angles
    between the diagonals of a parallelogram.
    """
    # Step 1: Define the vectors that generate the parallelogram.
    v1 = [6.0, 0.0, 0.0]
    v2 = [7.0/2.0, math.sqrt(13)/2.0, 0.0]

    # Step 2: Calculate the two diagonal vectors.
    d1 = [v1[i] + v2[i] for i in range(3)]
    d2 = [v1[i] - v2[i] for i in range(3)]

    # Step 3: Calculate the dot product of the diagonals.
    dot_product_diagonals = sum(d1[i] * d2[i] for i in range(3))

    # Step 4: Calculate the magnitudes of the diagonals.
    mag_d1 = math.sqrt(sum(x**2 for x in d1))
    mag_d2 = math.sqrt(sum(x**2 for x in d2))

    # Step 5: Calculate the cosine of the angle between the diagonal vectors.
    cos_angle_value = dot_product_diagonals / (mag_d1 * mag_d2)

    # Step 6: Determine all possible values for cos(theta).
    # The two angles at the intersection of the diagonals are supplementary.
    # Their cosines are equal in magnitude but opposite in sign.
    # We use a small tolerance for floating point comparison with zero.
    if abs(cos_angle_value) < 1e-9:
        possible_cos_values = [0.0]
    else:
        possible_cos_values = [cos_angle_value, -cos_angle_value]

    # Step 7: Calculate a and b.
    # a is the sum of all possible values of cos(theta).
    # b is the number of possible values of cos(theta).
    a = sum(possible_cos_values)
    b = len(possible_cos_values)

    # Step 8: Calculate the final result.
    result = a * b
    
    # Step 9: Print the steps of the final calculation as requested.
    val1 = possible_cos_values[0]
    # The 'final equation' is a * b = (sum of cos values) * (number of values) = result.
    # We show the numbers that form this equation.
    if b == 2:
        val2 = possible_cos_values[1]
        print(f"The possible values for cos(theta) are {val1:.4f} and {val2:.4f}.")
        print(f"a = sum of values = {val1:.4f} + ({val2:.4f}) = {a:.4f}")
        print(f"b = number of values = {b}")
        print(f"Final equation: ({val1:.4f} + {val2:.4f}) * {b} = {result:.3f}")
    else: # This case would apply if the diagonals were perpendicular.
        print(f"The only possible value for cos(theta) is {val1:.4f}.")
        print(f"a = sum of values = {a:.4f}")
        print(f"b = number of values = {b}")
        print(f"Final equation: {a:.4f} * {b} = {result:.3f}")


solve_parallelogram_diagonals()