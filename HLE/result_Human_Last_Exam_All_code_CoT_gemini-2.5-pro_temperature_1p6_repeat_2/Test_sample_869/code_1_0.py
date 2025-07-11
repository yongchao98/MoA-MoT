import math

def solve_probability():
    """
    Calculates the probability that for each individual, there exists a unique
    type of item for which they hold strictly more copies than any other individual.
    """

    # Step 1: Calculate S, the total number of ways to distribute the items.
    # There are 25 items in total, with 5 identical copies of 5 different types.
    # These are distributed among 5 individuals, 5 items each.
    # S = 25! / (5! * 5! * 5! * 5! * 5!) = 25! / (5!)^5
    try:
        s_value = math.factorial(25) // (math.factorial(5)**5)
    except OverflowError:
        print("Error: Calculation of S involves numbers too large to handle.")
        return

    # Step 2: Calculate F, the number of favorable distributions.
    # A favorable distribution requires each individual to have a unique "dominant" type.
    # The total number of ways to assign a unique dominant type to each individual is 5! (permutations).
    num_pairings = math.factorial(5)

    # We analyze configurations for a fixed pairing (e.g., individual 'j' is dominant in type 'j').
    # The total F is num_pairings times the sum of ways for these fixed-pairing configurations.
    # The configurations are broken into families based on the number of items 'd' on the dominant diagonal.
    # For a given configuration matrix C, the number of ways W(C) depends on the composition of each person's hand.

    # W_hand = 5! / (c1! * c2! * c3! * c4! * c5!), where ci is the count of items of type i.
    # W(C) = product over all individuals of their W_hand.

    # Case d=5: Hand is (5,0,0,0,0). Number of matrices is 1.
    # W_hand = 5!/5! = 1. W(C) = 1^5 = 1.
    ways_d5 = 1 * (1**5)

    # Case d=4: Hand is perm. of (4,1,0,0,0). Number of matrices is D_5=44.
    # W_hand = 5!/4! = 5. W(C) = 5^5.
    ways_d4 = 44 * (5**5)

    # Case d=3: Hand is perm. of (3,1,1,0,0). Number of matrices = 44.
    # W_hand = 5!/(3!*1!*1!) = 20. W(C) = 20^5.
    ways_d3 = 44 * (20**5)

    # Case d=2: Hand is perm. of (2,1,1,1,0). Number of matrices = 44.
    # W_hand = 5!/(2!*1!*1!*1!) = 60. W(C) = 60^5.
    ways_d2 = 44 * (60**5)
    
    # Sum of ways for a fixed pairing
    w_sum_fixed_pairing = ways_d5 + ways_d4 + ways_d3 + ways_d2
    
    # Total Favorable ways F
    f_value = num_pairings * w_sum_fixed_pairing

    # Step 3: Print the final equation and result.
    print("Calculation of the probability P = F / S")
    print("-" * 40)
    print(f"S (Total distributions) = 25! / (5!)^5 = {s_value}")
    print("")
    print("F (Favorable distributions) = 5! * (W_d5 + W_d4 + W_d3 + W_d2)")
    print(f"  where W_d5 (d=5 case) = 1 * (1^5) = {ways_d5}")
    print(f"        W_d4 (d=4 case) = 44 * (5^5) = {ways_d4}")
    print(f"        W_d3 (d=3 case) = 44 * (20^5) = {ways_d3}")
    print(f"        W_d2 (d=2 case) = 44 * (60^5) = {ways_d2}")
    print(f"Total F = {num_pairings} * ({ways_d5} + {ways_d4} + {ways_d3} + {ways_d2})")
    print(f"F = {f_value}")
    print("")
    print("Probability P = F / S")
    print(f"P = {f_value} / {s_value}")
    print(f"P = {f_value / s_value}")

solve_probability()