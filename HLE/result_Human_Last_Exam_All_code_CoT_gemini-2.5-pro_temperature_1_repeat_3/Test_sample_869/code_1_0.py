import math

def solve_probability():
    """
    Calculates the probability that each of five individuals receives a unique type of item
    in which they hold more copies than any other individual.
    """
    # N = 25 items, T = 5 types, 5 copies of each type.
    # Each of the 5 individuals gets 5 items.

    # S represents the total number of ways to distribute the items.
    # This is the number of unique permutations of the 25 items, which is given by
    # the multinomial coefficient 25! / (5! * 5! * 5! * 5! * 5!).
    try:
        fact_25 = math.factorial(25)
        fact_5 = math.factorial(5)
        S = fact_25 / (fact_5 ** 5)
    except OverflowError:
        print("Error: Calculation involves numbers too large to compute directly.")
        return

    # F represents the number of favorable distributions.
    # The constraints imply that the only favorable distributions are when each
    # individual receives 5 items of a single, unique type.
    # There are 5! ways to assign a unique type to each of the 5 individuals.
    # Each such assignment corresponds to exactly one permutation of the 25 items
    # (e.g., AAAAABBBBBCCCCCDDDDDEEEEE).
    # Therefore, the number of favorable outcomes F is 5!.
    F = math.factorial(5)

    # The probability P is F / S.
    P = F / S

    # Print the values for the equation P = F / S
    print("The total number of ways to distribute the items is:")
    print(f"S = 25! / (5!^5) = {S}")
    print("\nThe number of favorable distributions is:")
    print(f"F = 5! = {F}")
    print("\nThe probability P is F / S:")
    print(f"P = {F} / {S}")
    print(f"P = {P}")

solve_probability()
<<<120 * (120**5) / 2432902008176640000>>>