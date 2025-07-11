import math

def solve():
    """
    Solves the three parts of the user's question and prints the result
    in the specified format.
    """

    # Part (a): Calculation
    # The number is the sum of C(16,k) for k in {4, 8, 12, 16}.
    # The individual terms of the final equation are computed here.
    n_a = 16
    weights_a = [4, 8, 12, 16]
    terms_a = [math.comb(n_a, k) for k in weights_a]
    answer_a = sum(terms_a)

    # Part (b): Result from theoretical analysis
    answer_b = "yes"

    # Part (c): Calculation
    # The number is the sum of C(12,k) for k in {4, 8, 12}.
    # The individual terms of the final equation are computed here.
    n_c = 12
    weights_c = [4, 8, 12]
    terms_c = [math.comb(n_c, k) for k in weights_c]
    answer_c = sum(terms_c)

    # Print the final result in the requested format
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve()