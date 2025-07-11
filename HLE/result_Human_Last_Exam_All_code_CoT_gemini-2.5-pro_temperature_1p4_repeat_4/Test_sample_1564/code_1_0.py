import math

def solve_parallelogram_diagonals():
    """
    Calculates a*b based on the angles between the diagonals of a parallelogram.
    """
    # Step 1: Define the vectors u and v
    u = [6, 0, 0]
    v = [7/2, math.sqrt(13)/2, 0]

    # Step 2: Calculate the diagonal vectors d1 and d2
    d1 = [u[0] + v[0], u[1] + v[1], u[2] + v[2]]
    d2 = [u[0] - v[0], u[1] - v[1], u[2] - v[2]]

    # Step 3: Calculate the cosine of the angle between the diagonal vectors.
    # An efficient way to calculate the dot product d1 . d2 is using the identity:
    # (u+v).(u-v) = ||u||^2 - ||v||^2
    u_mag_sq = u[0]**2 + u[1]**2 + u[2]**2
    v_mag_sq = v[0]**2 + v[1]**2 + v[2]**2
    d1_dot_d2 = u_mag_sq - v_mag_sq

    # Calculate the magnitudes of the diagonals
    d1_mag = math.sqrt(d1[0]**2 + d1[1]**2 + d1[2]**2)
    d2_mag = math.sqrt(d2[0]**2 + d2[1]**2 + d2[2]**2)
    
    # Calculate one of the possible values for cos(theta)
    # This corresponds to the angle between the vectors d1 and d2.
    if d1_mag == 0 or d2_mag == 0:
      cos_alpha = 0.0
    else:
      cos_alpha = d1_dot_d2 / (d1_mag * d2_mag)

    # Step 4: Determine all possible values for cos(theta)
    # The two lines of the diagonals form two angles, alpha and (pi - alpha).
    # The cosines are cos(alpha) and -cos(alpha).
    # These are the possible values for cos(theta).
    # If cos_alpha is 0, the diagonals are perpendicular, and there is only one value (0).
    if cos_alpha == 0:
        possible_cos_values = [0.0]
    else:
        possible_cos_values = [cos_alpha, -cos_alpha]
        
    # Step 5: Calculate a and b
    # b is the number of possible values for cos(theta)
    b = len(possible_cos_values)
    # a is the sum of all possible values of cos(theta)
    a = sum(possible_cos_values)

    # Step 6: Calculate the final answer a * b and round it
    result = a * b
    
    print(f"The possible values for cos(theta) are: {possible_cos_values[0]:.4f} and {possible_cos_values[1]:.4f}")
    print(f"The number of possible values is b = {b}")
    print(f"The sum of all possible values is a = {a}")
    print(f"The final equation is a * b:")
    print(f"{a} * {b} = {result}")
    
    # Round the final result to the nearest thousandth
    rounded_result = round(result, 3)
    print(f"The result rounded to the nearest thousandth is {rounded_result:.3f}")

solve_parallelogram_diagonals()