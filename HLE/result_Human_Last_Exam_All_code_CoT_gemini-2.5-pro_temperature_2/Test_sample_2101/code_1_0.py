import math

def solve_escape_probability():
    """
    Calculates the difference in escape probabilities for an electron in an isosceles right triangle.

    The solution is based on the principle that the escape probability through a given side
    of a triangle is proportional to the angle at the opposite vertex.
    """

    # 1. Define the internal angles of the isosceles right triangle in radians.
    # The angle opposite the hypotenuse is 90 degrees.
    angle_opp_hypotenuse = math.pi / 2
    # The angles opposite the two equal legs are 45 degrees each.
    angle_opp_leg = math.pi / 4
    
    # The sum of all angles is the constant of proportionality for the probabilities.
    total_angle_sum = angle_opp_hypotenuse + angle_opp_leg + angle_opp_leg

    # 2. Calculate the escape probabilities.
    # P(side) = (angle opposite side) / (sum of all angles)
    prob_hypotenuse = angle_opp_hypotenuse / total_angle_sum
    prob_one_leg = angle_opp_leg / total_angle_sum
    prob_two_legs = prob_one_leg + prob_one_leg
    
    # 3. Calculate the difference.
    difference = prob_hypotenuse - prob_two_legs
    
    # 4. Print the results step-by-step, including the final equation.
    print("This problem can be solved using a principle of geometric probability.")
    print("The probability of escaping through a side is proportional to the angle at the opposite vertex.\n")
    print(f"Angle opposite hypotenuse: 90 degrees (pi/2 radians)")
    print(f"Angle opposite each leg: 45 degrees (pi/4 radians)\n")
    
    print(f"Probability of escaping through the hypotenuse = (pi/2) / (pi/2 + pi/4 + pi/4) = {prob_hypotenuse:.1f}")
    print(f"Probability of escaping through one leg = (pi/4) / (pi/2 + pi/4 + pi/4) = {prob_one_leg:.2f}")
    print(f"Probability of escaping through either leg = {prob_one_leg:.2f} + {prob_one_leg:.2f} = {prob_two_legs:.1f}\n")
    
    print("The difference is the probability of escaping through the hypotenuse minus the probability of escaping through the legs.")
    print(f"Difference = P(hypotenuse) - (P(leg1) + P(leg2))")
    print(f"Difference = {prob_hypotenuse:.1f} - ({prob_one_leg:.2f} + {prob_one_leg:.2f})")
    print(f"Difference = {prob_hypotenuse:.1f} - {prob_two_legs:.1f} = {difference:.1f}")

solve_escape_probability()