import math

def solve_electron_escape_problem():
    """
    Calculates the difference between the probability that an electron escapes
    through the hypotenuse and the probability that it escapes through either
    of the two legs of an isosceles right triangle.

    The solution is based on a theorem from integral geometry for the average
    angle subtended by a side of a triangle from a random interior point.
    """

    # The internal angles of an isosceles right triangle
    # Let the vertices be A, B, C, with the right angle at B.
    angle_A = math.pi / 4  # 45 degrees
    angle_B = math.pi / 2  # 90 degrees
    angle_C = math.pi / 4  # 45 degrees

    # The average angle subtended by a side is the sum of the internal angles at its endpoints.
    # Average angle for the hypotenuse (side AC)
    avg_angle_hypotenuse = angle_A + angle_C

    # Average angle for the first leg (side AB)
    avg_angle_leg1 = angle_A + angle_B

    # Average angle for the second leg (side BC)
    avg_angle_leg2 = angle_B + angle_C
    
    # The probability of escaping through a side is the average subtended angle divided by 2*pi.
    # Probability of escaping through the hypotenuse
    prob_hypotenuse = avg_angle_hypotenuse / (2 * math.pi)

    # Probability of escaping through the legs
    prob_leg1 = avg_angle_leg1 / (2 * math.pi)
    prob_leg2 = avg_angle_leg2 / (2 * math.pi)
    prob_legs = prob_leg1 + prob_leg2

    # The difference between the probabilities
    difference = prob_hypotenuse - prob_legs

    print("Step 1: Calculate the probability of escaping through the hypotenuse (P_H).")
    print(f"P_H = (Angle A + Angle C) / (2 * pi) = ({angle_A/math.pi:.2f}*pi + {angle_C/math.pi:.2f}*pi) / (2*pi) = {prob_hypotenuse}")
    print("\nStep 2: Calculate the probability of escaping through the legs (P_L).")
    print(f"P_L = P_L1 + P_L2")
    print(f"P_L1 = (Angle A + Angle B) / (2 * pi) = ({angle_A/math.pi:.2f}*pi + {angle_B/math.pi:.2f}*pi) / (2*pi) = {prob_leg1}")
    print(f"P_L2 = (Angle B + Angle C) / (2 * pi) = ({angle_B/math.pi:.2f}*pi + {angle_C/math.pi:.2f}*pi) / (2*pi) = {prob_leg2}")
    print(f"P_L = {prob_leg1} + {prob_leg2} = {prob_legs}")
    print("\nStep 3: Calculate the difference.")
    print(f"Difference = P_H - P_L")
    print(f"Difference = {prob_hypotenuse} - {prob_legs} = {difference}")

solve_electron_escape_problem()
<<< -0.5 >>>