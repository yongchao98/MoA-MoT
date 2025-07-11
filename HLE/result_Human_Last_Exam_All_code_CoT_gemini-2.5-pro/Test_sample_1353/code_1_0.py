import math

def solve_dh_problem():
    """
    Solves the three-part problem about diagonal harmonics.
    """

    # Part a
    # A string starter P has bi-degree (a, b) = (4, 3).
    # The sl(2) string is symmetric. The starting (lowest weight) bi-degree is (a, b).
    # The terminal (highest weight) bi-degree is (b, a).
    a_val, b_val = 4, 3
    terminal_bidegree = (b_val, a_val)
    answer_a = f"({terminal_bidegree[0]}, {terminal_bidegree[1]})"

    # Part b
    # For a polynomial of bi-degree (a, b) to be a string starter constructed
    # with b indices r_1, ..., r_b, these indices must form a partition of b.
    # The condition r_1 + ... + r_b = b for positive integers r_i implies r_i=1 for all i.
    # The x-degree 'a' is then given by the statistic n(lambda) for lambda=(1,1,...,1).
    # a = sum_{i=1 to b} (i-1)*r_i = sum_{i=1 to b} (i-1) = b*(b-1)/2
    # So the condition is a = C(b, 2).
    # We will output the equation a = b(b-1)/2, showing the numbers.
    b_sym = 'b'
    a_sym = 'a'
    # Using C(b, 2) is also a valid expression.
    # a = b * (b - 1) / 2
    answer_b = f"a = {b_sym}({b_sym}-1)/2 = (({b_sym}^2 - {b_sym})/2)"

    # Part c
    # Can a polynomial of bi-degree (5, 2) be constructed using E_{r,0}?
    # Interpreting E_{r,0} as multiplication by p_r(Y) = sum(y_i^r).
    # A polynomial constructed this way, e.g., p_r(Y) * Q, is not in DH_n
    # because it is not annihilated by the defining differential operators.
    # For example, sum(d/dy_i) applied to p_1(Y) gives n, not 0.
    # Therefore, it is not possible.
    answer_c = "No"

    print(f"a) {answer_a}")
    print(f"b) {answer_b}")
    print(f"c) {answer_c}")

solve_dh_problem()
<<<a) (3, 4) b) a = b(b-1)/2 = ((b^2 - b)/2) c) No>>>