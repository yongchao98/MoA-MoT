import math

def solve_electron_escape_probability():
    """
    Calculates the difference in escape probabilities for an electron in an isosceles right triangle.
    """
    pi = math.pi

    # Step 1: Define the internal angles of the isosceles right triangle
    angle_opposite_hypotenuse = pi / 2  # 90 degrees
    angle_opposite_leg = pi / 4        # 45 degrees

    print("Step 1: The internal angles of the triangle are 90, 45, and 45 degrees.")
    print(f"In radians, the angle opposite the hypotenuse is {angle_opposite_hypotenuse/pi:.2f}π")
    print(f"In radians, the angle opposite each leg is {angle_opposite_leg/pi:.2f}π\n")

    # Step 2: Use the theorem <A_s> = π - α to find the average subtended angles.
    # <A_s> is the average angle subtended by a side s.
    # α is the triangle's internal angle at the vertex opposite side s.
    avg_angle_hypotenuse = pi - angle_opposite_hypotenuse
    avg_angle_leg = pi - angle_opposite_leg
    
    print("Step 2: Calculate the average angle subtended by each side from a random point inside.")
    print("We use the theorem: Average Subtended Angle = π - (Opposite Internal Angle)")
    print(f"Average angle for hypotenuse = π - π/2 = {avg_angle_hypotenuse/pi:.2f}π")
    print(f"Average angle for a leg = π - π/4 = {avg_angle_leg/pi:.2f}π\n")

    # Step 3: Calculate the probability of escape through each side.
    # P(s) = <A_s> / 2π
    prob_hypotenuse = avg_angle_hypotenuse / (2 * pi)
    prob_leg = avg_angle_leg / (2 * pi)
    
    print("Step 3: Calculate the escape probability for each side.")
    print("We use the formula: Probability = Average Subtended Angle / 2π")
    print(f"Probability of escaping through the hypotenuse = ({avg_angle_hypotenuse/pi:.2f}π) / 2π = {prob_hypotenuse}")
    print(f"Probability of escaping through one leg = ({avg_angle_leg/pi:.2f}π) / 2π = {prob_leg}\n")

    # Step 4: Calculate the total probability for the legs and find the final difference.
    prob_both_legs = 2 * prob_leg
    difference = prob_hypotenuse - prob_both_legs

    print("Step 4: Calculate the final difference.")
    print("The quantity we want is P(hypotenuse) - P(either leg).")
    print(f"P(either leg) = P(leg 1) + P(leg 2) = {prob_leg} + {prob_leg} = {prob_both_legs}")
    print(f"Difference = P(hypotenuse) - P(either leg)")
    print(f"Difference = {prob_hypotenuse} - {prob_both_legs} = {difference}")

solve_electron_escape_probability()