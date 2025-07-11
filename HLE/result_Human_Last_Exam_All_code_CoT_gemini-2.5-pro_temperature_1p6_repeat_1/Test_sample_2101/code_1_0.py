import math

def solve_electron_escape_probability():
    """
    Calculates the difference between the probability that an electron escapes
    through the hypotenuse and the probability that it escapes through either
    of the two legs of an isosceles right triangle.

    The solution is based on the principle that for a particle generated at a
    uniform random point and moving in a random direction inside a convex
    2D shape, the probability of escaping through any part of the boundary
    is proportional to the length of that part.
    """

    # We can set the leg length L to any value, as it will cancel out. Let's use 1.
    L = 1.0

    # Calculate the lengths of the sides
    leg_length_total = 2 * L
    hypotenuse_length = L * math.sqrt(2)

    # Calculate the total perimeter
    perimeter = leg_length_total + hypotenuse_length

    # Calculate the probabilities
    prob_hypotenuse = hypotenuse_length / perimeter
    prob_legs = leg_length_total / perimeter

    # Calculate the difference
    difference = prob_hypotenuse - prob_legs
    
    # The final equation is P(hypotenuse) - P(legs) = difference.
    # Let's print the equation with the calculated numbers.
    print("The probability of escape is proportional to the length of the boundary.")
    print(f"Let the leg length be {L}.")
    print(f"The total length of the two legs is {leg_length_total:.4f}.")
    print(f"The length of the hypotenuse is sqrt(2) * {L} = {hypotenuse_length:.4f}.")
    print(f"The total perimeter is {leg_length_total:.4f} + {hypotenuse_length:.4f} = {perimeter:.4f}.")
    print("-" * 30)
    print("The equation for the difference in probabilities is:")
    print("P(hypotenuse) - P(legs)")
    print(f"= ({hypotenuse_length:.4f} / {perimeter:.4f}) - ({leg_length_total:.4f} / {perimeter:.4f})")
    print(f"= {prob_hypotenuse:.4f} - {prob_legs:.4f}")
    
    # The simplified exact form is (sqrt(2) - 2) / (sqrt(2) + 2), which is 2*sqrt(2) - 3.
    exact_difference = 2 * math.sqrt(2) - 3
    print(f"= {exact_difference:.4f}")
    
    # The final expression we are calculating:
    print("\nFinal calculation:")
    print(f"{prob_hypotenuse:.4f} - {prob_legs:.4f} = {difference:.4f}")

if __name__ == '__main__':
    solve_electron_escape_probability()
