import math

def h(x, y):
    """
    The potential kernel h(z) for the 2D random walk is asymptotically
    proportional to log|z|. We use this approximation.
    The state space is Z^2 \ {0}.
    """
    dist_sq = float(x**2 + y**2)
    if dist_sq == 0:
        # This point is not in the state space of the conditioned walk.
        return -math.inf
    return math.log(math.sqrt(dist_sq))

def demonstrate_weakening_bias():
    """
    This function calculates and prints the ratio of outward to inward
    transition probabilities for the h-transformed random walk at various
    distances from the origin.
    """
    # A list of distances (radii) from the origin to test.
    radii = [10, 100, 1000, 10000, 100000]

    print("This script demonstrates that the outward drift of the conditioned random walk")
    print("weakens as the walk moves away from the origin. This supports the argument")
    print("that an infinite but very sparse set can be transient.\n")
    print("We calculate the following equation for various R:")
    print("  Bias Ratio = P(move out) / P(move in)\n")

    print(f"{'Distance R':>12} | {'Bias Ratio':>20}")
    print("-" * 35)

    for R in radii:
        if R <= 1:
            continue

        # Consider a point x = (R, 0)
        # Its outward neighbor is y_out = (R+1, 0)
        # Its inward neighbor is y_in = (R-1, 0)

        # The transition probabilities for the h-transformed walk P_hat are related by:
        # P_hat(x, y_out) / P_hat(x, y_in) = h(y_out) / h(y_in)
        # since the simple random walk probabilities P(x, y_out) and P(x, y_in) are equal.
        
        h_yout = h(R + 1, 0)
        h_yin = h(R - 1, 0)
        
        # This is the final equation we compute:
        if h_yin == 0: # Avoid division by zero if R=1
            continue
        bias_ratio = h_yout / h_yin

        # We output each number R and the resulting bias_ratio from our equation.
        print(f"{R:12d} | {bias_ratio:20.15f}")

    print("\nAs R increases, the ratio approaches 1. This means the walk behaves")
    print("more like an unbiased simple random walk far from the origin.")

if __name__ == '__main__':
    demonstrate_weakening_bias()
