from fractions import Fraction

def solve_coin_probability():
    """
    Calculates the probability of getting all heads on three coins,
    given that the number of tails is always even.
    """

    # 1. Define initial probabilities as fractions for precision
    p_h = Fraction(1, 3)
    p_t = Fraction(2, 3)

    print("Step 1: Define the base probabilities.")
    print(f"The probability of a single coin turning heads, P(H) = {p_h}")
    print(f"The probability of a single coin turning tails, P(T) = {p_t}")
    print("-" * 30)

    # 2. Calculate the probability of the event we are interested in: all heads (HHH).
    # This is also the case for "0 tails", which is an even number.
    prob_all_heads = p_h ** 3
    print("Step 2: Calculate the probability of getting 0 tails (all heads: HHH).")
    print(f"P(0 Tails) = P(H) * P(H) * P(H) = ({p_h}) * ({p_h}) * ({p_h}) = {prob_all_heads}")
    print("-" * 30)

    # 3. Calculate the probability of the other case satisfying the condition: 2 tails.
    # There are 3 combinations for 2 tails and 1 head (HTT, THT, TTH).
    num_combinations_2_tails = 3
    prob_one_combination_2_tails = p_h * (p_t ** 2)
    prob_total_2_tails = num_combinations_2_tails * prob_one_combination_2_tails
    print("Step 3: Calculate the probability of getting 2 tails.")
    print("This can happen in 3 ways (HTT, THT, TTH).")
    print(f"The probability of one specific combination (e.g., HTT) is P(H)*P(T)*P(T) = {prob_one_combination_2_tails}")
    print(f"The total probability for 2 tails is 3 * {prob_one_combination_2_tails} = {prob_total_2_tails}")
    print("-" * 30)

    # 4. Calculate the total probability of the condition (even number of tails).
    # This is P(0 Tails) + P(2 Tails).
    prob_even_tails = prob_all_heads + prob_total_2_tails
    print("Step 4: Calculate the total probability of the number of tails being even.")
    print(f"P(Even Tails) = P(0 Tails) + P(2 Tails) = {prob_all_heads} + {prob_total_2_tails} = {prob_even_tails}")
    print("-" * 30)

    # 5. Calculate the final conditional probability.
    # P(All Heads | Even Tails) = P(All Heads) / P(Even Tails)
    # Note: P(All Heads and Even Tails) is just P(All Heads) because
    # the 'All Heads' outcome (0 tails) already satisfies the 'Even Tails' condition.
    final_probability = prob_all_heads / prob_even_tails
    print("Step 5: Calculate the final conditional probability.")
    print("P(All Heads | Even Tails) = P(All Heads) / P(Even Tails)")
    print("\nThe final equation is:")
    print(f"{prob_all_heads} / {prob_even_tails} = {final_probability}")

solve_coin_probability()
<<<1/13>>>