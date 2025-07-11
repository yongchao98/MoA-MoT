import math

def solve_coin_problem():
    """
    Calculates the conditional probability of getting all heads given that the number of tails is even for 3 biased coins.
    """
    # Step 1: Define probabilities for a single coin
    p_heads = 1/3
    p_tails = 2/3

    # Step 2: Calculate the probability of the condition (even number of tails)
    # This can be 0 tails (HHH) or 2 tails (HTT, THT, TTH).

    # Probability of 0 tails (HHH)
    prob_0_tails = p_heads ** 3
    num_0_tails = int(1/prob_0_tails)

    # Probability of 2 tails (e.g., HTT)
    # There are 3 ways to get 2 tails: C(3,2) = 3
    num_ways_2_tails = 3
    prob_2_tails = num_ways_2_tails * (p_heads**1) * (p_tails**2)
    # Using fractions for precision in the printout
    numerator_2_tails = num_ways_2_tails * (2**2)
    denominator_2_tails = 3**3

    # Total probability of the condition (even tails)
    prob_even_tails = prob_0_tails + prob_2_tails
    # Using fractions for precision in the printout
    numerator_even_tails = 1 + numerator_2_tails
    denominator_even_tails = 3**3


    # Step 3: Identify the probability of the event we're interested in AND the condition being true.
    # Event: All heads (HHH).
    # This event (0 tails) satisfies the condition of having an even number of tails.
    # So, P(All Heads and Even Tails) is just P(All Heads).
    prob_all_heads_and_even_tails = prob_0_tails
    num_all_heads_and_even_tails = 1 # numerator for 1/27

    # Step 4: Apply the conditional probability formula
    # P(All Heads | Even Tails) = P(All Heads and Even Tails) / P(Even Tails)
    conditional_prob = prob_all_heads_and_even_tails / prob_even_tails

    # Step 5: Print the explanation and the final equation
    print("Let P(H) be the probability of heads and P(T) be the probability of tails.")
    print(f"P(H) = 1/3, P(T) = 2/3\n")
    print("The condition is that the number of tails is even (0 or 2 tails).")
    print(f"P(0 Tails) = P(HHH) = (1/3)^3 = 1/{num_0_tails}")
    print(f"P(2 Tails) = 3 * P(H) * P(T)^2 = 3 * (1/3) * (2/3)^2 = {numerator_2_tails}/{denominator_2_tails}")
    print(f"P(Even Tails) = P(0 Tails) + P(2 Tails) = 1/{num_0_tails} + {numerator_2_tails}/{denominator_2_tails} = {numerator_even_tails}/{denominator_even_tails}\n")
    print("The event of interest is getting all heads (HHH).")
    print(f"P(All Heads) = 1/{num_0_tails}\n")
    print("The probability of 'All Heads' given 'Even Tails' is:")
    print(f"P(All Heads | Even Tails) = P(All Heads) / P(Even Tails)")
    print(f"= (1/{num_0_tails}) / ({numerator_even_tails}/{denominator_even_tails})")
    print(f"= {int(conditional_prob**-1)}^-1 = {conditional_prob:.4f}\n")
    print("Final Equation:")
    print(f"1/{num_all_heads_and_even_tails} / {numerator_even_tails} = {int(conditional_prob**-1)}")
    print(f"The final probability is 1/{int(conditional_prob**-1)}")


solve_coin_problem()
<<<1/13>>>