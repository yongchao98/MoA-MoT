import sys

def solve_coin_problem():
    """
    Calculates the probability of getting all heads on three coin tosses,
    given that the number of tails is even.
    P(H) = 1/3.
    """
    # Numerators and denominators to work with fractions
    p_h_num = 1
    p_h_den = 3
    p_t_num = 2
    p_t_den = 3
    
    # Event A: All three coins are heads (HHH)
    # P(A) = (1/3)^3 = 1/27
    p_A_num = p_h_num ** 3
    p_A_den = p_h_den ** 3

    # Event B: The number of tails is even (0 or 2 tails)
    # This is the sum of P(0 tails) and P(2 tails)

    # Probability of 0 tails (HHH), which is the same as P(A)
    p_0_tails_num = p_A_num
    p_0_tails_den = p_A_den

    # Probability of 2 tails (e.g., HTT)
    # P(HTT) = (1/3)*(2/3)*(2/3) = 4/27
    # There are 3 combinations for 2 tails (HTT, THT, TTH)
    # P(2 tails) = 3 * P(HTT) = 3 * (4/27) = 12/27
    p_2_tails_num = 3 * p_h_num * (p_t_num ** 2)
    p_2_tails_den = p_h_den * (p_t_den ** 2)

    # Probability of Event B = P(0 tails) + P(2 tails)
    # P(B) = 1/27 + 12/27 = 13/27
    p_B_num = p_0_tails_num + p_2_tails_num
    p_B_den = p_A_den # Common denominator

    # We want to find P(A | B) = P(A and B) / P(B)
    # The event "A and B" (all heads and an even number of tails) is just event A.
    # So, P(A and B) = P(A).
    # P(A|B) = P(A) / P(B)
    final_numerator = p_A_num
    final_denominator = p_B_num

    print("Let A be the event 'all heads' and B be the event 'even number of tails'.")
    print("We want to calculate P(A|B) = P(A and B) / P(B).\n")

    print(f"The probability of event A (all heads) is P(A) = {p_A_num}/{p_A_den}.")
    print("The event B (even tails) means 0 tails or 2 tails.")
    print(f"P(0 tails) = P(HHH) = {p_0_tails_num}/{p_0_tails_den}")
    print(f"P(2 tails) = P(HTT, THT, TTH) = {p_2_tails_num}/{p_2_tails_den}")
    print(f"So, P(B) = P(0 tails) + P(2 tails) = {p_0_tails_num}/{p_0_tails_den} + {p_2_tails_num}/{p_2_tails_den} = {p_B_num}/{p_B_den}.\n")
    
    print("Since 'A and B' is the same as A, P(A|B) = P(A) / P(B).\n")

    print("The final equation with the calculated probabilities is:")
    print(f"({p_A_num}/{p_A_den}) / ({p_B_num}/{p_B_den})\n")
    
    print(f"This simplifies to {final_numerator}/{final_denominator}.")

    # Use 'write' to avoid extra newline for the final answer format
    sys.stdout.write(f"<<<{final_numerator}/{final_denominator}>>>")

solve_coin_problem()