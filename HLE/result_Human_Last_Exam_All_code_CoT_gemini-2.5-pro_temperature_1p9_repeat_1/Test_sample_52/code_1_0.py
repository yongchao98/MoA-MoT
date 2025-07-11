from fractions import Fraction

def solve_coin_problem():
    """
    Calculates the probability of three coins all turning heads,
    given that the number of tails is even.
    """
    # Step 1: Define the base probabilities for a single coin.
    p_h = Fraction(1, 3)
    p_t = 1 - p_h

    print(f"The probability of a single coin turning Heads is P(H) = {p_h}.")
    print(f"The probability of a single coin turning Tails is P(T) = {p_t}.")
    print("-" * 30)

    # Step 2: Calculate the probability of the condition (B: even number of tails).
    # This means 0 tails (HHH) or 2 tails (HTT, THT, TTH).

    # Probability of 0 tails (HHH)
    p_0_tails = p_h ** 3

    # Probability of 2 tails. There are 3 combinations (HTT, THT, TTH).
    # Each combination has a probability of p_h * p_t**2.
    p_2_tails = 3 * p_h * (p_t ** 2)

    # Total probability of the condition (P(B))
    p_even_tails = p_0_tails + p_2_tails

    print("The condition is that the number of tails is even (0 or 2 tails).")
    print(f"Probability of 0 tails (HHH): {p_h}*{p_h}*{p_h} = {p_0_tails}")
    print(f"Probability of 2 tails (any order): 3 * {p_h} * {p_t}**2 = {p_2_tails}")
    print(f"Total probability of an even number of tails is P(B) = {p_0_tails} + {p_2_tails} = {p_even_tails}")
    print("-" * 30)

    # Step 3: Calculate the probability of the intersection (A and B).
    # Event A is "all heads". The intersection of "all heads" and "even number of tails"
    # is just "all heads", as 0 tails is an even number.
    p_all_heads_and_even_tails = p_0_tails

    print("We want to find the probability of 'all heads' GIVEN 'even number of tails'.")
    print("Using the formula P(A|B) = P(A and B) / P(B).")
    print(f"P(A and B) is the probability of 'all heads', which is {p_all_heads_and_even_tails}.")
    print("-" * 30)

    # Step 4: Calculate the final conditional probability.
    final_prob = p_all_heads_and_even_tails / p_even_tails

    print("The final calculation is:")
    # Using format to explicitly show each number in the equation.
    print(f"({p_all_heads_and_even_tails}) / ({p_even_tails}) = {final_prob}")

if __name__ == '__main__':
    solve_coin_problem()
<<<1/13>>>