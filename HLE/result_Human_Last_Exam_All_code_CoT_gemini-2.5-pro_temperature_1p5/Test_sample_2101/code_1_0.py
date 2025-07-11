import math

def solve_escape_probability():
    """
    Calculates the difference in escape probabilities for an electron from an
    isosceles right triangle.
    """
    # Step 1 & 2: Define geometry and calculate side lengths.
    # We can set the leg length L=1.0 without loss of generality, as the L factor
    # will cancel out when calculating probability ratios.
    leg_length = 1.0
    
    # The triangle has two legs of equal length and one hypotenuse.
    len_leg_1 = leg_length
    len_leg_2 = leg_length
    len_hypotenuse = math.sqrt(leg_length**2 + leg_length**2) # This is leg_length * sqrt(2)

    # Step 3: Calculate the total perimeter.
    perimeter = len_leg_1 + len_leg_2 + len_hypotenuse

    # Step 4: Apply the principle that escape probability is proportional to side length.
    # P(escape through a side) = (length of side) / (total perimeter)
    # This holds for any convex shape given the problem's assumptions.

    # Step 5: Calculate the probability of escaping through the hypotenuse.
    prob_hypotenuse = len_hypotenuse / perimeter

    # Step 6: Calculate the probability of escaping through either of the two legs.
    prob_legs = (len_leg_1 + len_leg_2) / perimeter

    # Step 7: Compute the final difference.
    difference = prob_hypotenuse - prob_legs

    # Output the results and the final equation as requested.
    print("Based on principles of geometric probability, the escape chance through any side is its length divided by the total perimeter.")
    print("-" * 20)
    print(f"Assuming leg lengths are {leg_length:.1f}:")
    print(f"Hypotenuse length: sqrt({leg_length:.1f}^2 + {leg_length:.1f}^2) = {len_hypotenuse:.4f}")
    print(f"Total Perimeter: {len_leg_1:.1f} + {len_leg_2:.1f} + {len_hypotenuse:.4f} = {perimeter:.4f}")
    print("-" * 20)
    
    print(f"The probability of escaping through the hypotenuse is {len_hypotenuse:.4f} / {perimeter:.4f} = {prob_hypotenuse:.4f}")
    print(f"The probability of escaping through the legs is ({len_leg_1:.1f} + {len_leg_2:.1f}) / {perimeter:.4f} = {prob_legs:.4f}")
    print("-" * 20)
    
    print("The final equation for the difference is:")
    print(f"P(hypotenuse) - P(legs) = {prob_hypotenuse:.4f} - {prob_legs:.4f}")
    print(f"Result = {difference:.4f}")
    
    # The precise analytical result is 2*sqrt(2) - 3
    print(f"\nThe exact value of this difference is 2*sqrt(2) - 3, which is {2*math.sqrt(2)-3}.")

solve_escape_probability()
