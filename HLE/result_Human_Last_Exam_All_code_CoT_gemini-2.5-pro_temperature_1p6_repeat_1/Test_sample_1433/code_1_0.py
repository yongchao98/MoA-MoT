import math

def solve():
    """
    Solves the user's questions about random walks on a hypercube.
    """

    # Part 1: Expected time to meet on a 14-dimensional hypercube (EX_14)
    d_14 = 14
    # For an even dimension d, the expected meeting time is 2^(d-1).
    ex_14 = 2**(d_14 - 1)
    print(f"The expected time they meet on a 14-hypercube is EX_14 = {int(ex_14)}")

    # Part 2: Variance of the meeting time (D^2 X_14)
    # The variance is approximately the square of the expectation for large d.
    var_x_14 = ex_14**2
    print(f"The variance of the time they meet is D^2X_14 = {int(var_x_14)}")

    # Part 3: Expected time to meet on a 15-dimensional hypercube (EX_15)
    # For an odd dimension d, the parity of the distance between the walkers is always odd.
    # Since the starting distance is d (odd), it can never become 0. They never meet.
    print(f"The expected time they meet on a 15-hypercube is EX_15 = inf")

    # Part 4: Check the inequality for even d
    # Is it true that EX_d <= (d/2) * (d^d / d!) for even d?
    # Let's test for d = 14.
    lhs = ex_14
    rhs = (d_14 / 2) * (d_14**d_14 / math.factorial(d_14))
    is_true = "yes" if lhs <= rhs else "no"
    print(f"Is it true that for even d, EX_d <= (d/2) * (d^d / d!)? {is_true}")
    
solve()
<<<8192, 67108864, inf, yes>>>