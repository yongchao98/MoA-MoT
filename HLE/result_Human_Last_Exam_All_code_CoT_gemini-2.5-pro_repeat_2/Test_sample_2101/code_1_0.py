import math

def solve_escape_probability():
    """
    Calculates the difference in escape probabilities for an electron in an
    isosceles right triangle.
    """
    pi = math.pi

    # Step 1: Define the angles of the isosceles right triangle in radians.
    # The hypotenuse is opposite the 90-degree angle.
    # The two legs are opposite the 45-degree angles.
    angle_opposite_hypotenuse = pi / 2.0  # 90 degrees
    angle_opposite_leg = pi / 4.0       # 45 degrees

    print("Step 1: Define the triangle's angles (in radians)")
    print(f"Angle opposite the hypotenuse = {angle_opposite_hypotenuse:.4f} (pi/2)")
    print(f"Angle opposite a leg = {angle_opposite_leg:.4f} (pi/4)")
    print("-" * 30)

    # Step 2: Use the formula P(side) = (pi - alpha) / (2 * pi)
    # where alpha is the angle opposite the side.
    print("Step 2: Calculate the probability of escaping through the hypotenuse (P_hyp)")
    # P_hyp = (pi - angle_opposite_hypotenuse) / (2 * pi)
    # P_hyp = (pi - pi/2) / (2*pi) = (pi/2) / (2*pi) = 1/4
    p_hyp = (pi - angle_opposite_hypotenuse) / (2 * pi)
    print(f"P_hyp = (pi - {angle_opposite_hypotenuse:.4f}) / (2 * pi) = {p_hyp}")
    print("-" * 30)

    print("Step 3: Calculate the probability of escaping through one leg (P_leg)")
    # P_leg = (pi - angle_opposite_leg) / (2 * pi)
    # P_leg = (pi - pi/4) / (2*pi) = (3*pi/4) / (2*pi) = 3/8
    p_leg = (pi - angle_opposite_leg) / (2 * pi)
    print(f"P_leg = (pi - {angle_opposite_leg:.4f}) / (2 * pi) = {p_leg}")
    print("-" * 30)
    
    # Step 4: The probability of escaping through either leg is 2 * P_leg
    print("Step 4: Calculate the total probability of escaping through the legs (P_legs)")
    p_legs = 2 * p_leg
    print(f"P_legs = P_leg + P_leg = {p_leg} + {p_leg} = {p_legs}")
    print("-" * 30)

    # Step 5: Find the final difference
    print("Step 5: Calculate the final difference: P_hyp - P_legs")
    difference = p_hyp - p_legs
    print(f"Difference = P_hyp - P_legs = {p_hyp} - {p_legs} = {difference}")
    print("-" * 30)
    
    # Final answer
    print("The final result is -0.5 or -1/2.")

solve_escape_probability()
