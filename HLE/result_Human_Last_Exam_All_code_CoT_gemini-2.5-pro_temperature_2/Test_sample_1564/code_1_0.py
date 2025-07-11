import math

def main():
    """
    Calculates the product of the sum and number of possible cosines of the angle between a parallelogram's diagonals.
    """

    # 1. Define the generating vectors
    u = [6.0, 0.0, 0.0]
    v = [7.0/2.0, math.sqrt(13)/2.0, 0.0]

    # 2. Calculate the diagonal vectors
    d1 = [u[0] + v[0], u[1] + v[1], u[2] + v[2]]
    d2 = [u[0] - v[0], u[1] - v[1], u[2] - v[2]]

    # Helper function for dot product
    def dot_product(v1, v2):
        return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

    # Helper function for magnitude
    def magnitude(v):
        return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)

    # 3. Analyze possible values for cos(theta)
    # The dot product of the diagonals is needed.
    d1_dot_d2 = dot_product(d1, d2)
    
    # Check for orthogonality (if dot product is zero)
    # Use a small tolerance for floating point comparison
    if abs(d1_dot_d2) < 1e-9:
        possible_cos_values = [0.0]
    else:
        # 4. Calculate one possible cosine value, c
        mag_d1 = magnitude(d1)
        mag_d2 = magnitude(d2)
        c = d1_dot_d2 / (mag_d1 * mag_d2)
        possible_cos_values = [c, -c]

    # 5. Determine 'a' (sum of values) and 'b' (number of values)
    a = sum(possible_cos_values)
    b = len(possible_cos_values)

    # 6. Compute the final result
    result = a * b

    # Round the final answer to the nearest thousandth
    final_answer = round(result, 3)

    # Output the steps and the final equation
    print("The problem is to find a * b, where 'a' is the sum of all possible values of cos(theta) and 'b' is the number of possible values.")
    print("The two possible values for cos(theta) are c and -c (unless the diagonals are perpendicular).")
    print(f"The list of possible values for cos(theta) is: {possible_cos_values}")
    
    print(f"\na, the sum of these values, is: {a:.3f}")
    print(f"b, the number of these values, is: {b}")

    # Output each number in the final equation
    print(f"\nThe final calculation is a * b.")
    print(f"{a:.3f} * {b} = {result:.3f}")

    print(f"\nThe answer rounded to the nearest thousandth is {final_answer:.3f}.")

    print(f"\n<<<{final_answer:.3f}>>>")

if __name__ == '__main__':
    main()