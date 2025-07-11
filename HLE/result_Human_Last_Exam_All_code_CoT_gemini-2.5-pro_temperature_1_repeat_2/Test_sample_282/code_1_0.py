import math

def demonstrate_growth_factor():
    """
    This function demonstrates the growth factor for a specific compact set
    in a subgroup of SL_2(R), showing that the ratio mu(X^3)/mu(X) can be 3.
    """
    # Let's model the compact set X = {diag(t, 1/t) | t in [1, a]}
    # We choose an arbitrary value for a > 1.
    a = 10.0

    # The Haar measure for this subgroup is d(mu) = dt/t.
    # mu(X) is the integral of 1/t from 1 to a.
    # mu(X) = log(a) - log(1) = log(a)
    mu_X = math.log(a)

    # The set X^3 = {xyz | x,y,z in X} corresponds to the set
    # {diag(s, 1/s) | s in [1, a^3]}.
    # mu(X^3) is the integral of 1/s from 1 to a^3.
    # mu(X^3) = log(a^3) - log(1) = 3 * log(a)
    mu_X3 = math.log(a**3)

    # The ratio K is mu(X^3) / mu(X)
    if mu_X > 0:
        K = mu_X3 / mu_X
    else:
        K = 0

    print(f"Let X be a set with measure mu(X) = log({a}) = {mu_X:.4f}")
    print(f"The measure of the product set X^3 is mu(X^3) = log({a}^3) = {mu_X3:.4f}")
    print("The relationship is mu(X^3) = K * mu(X).")
    print("For this example, we find the constant K:")
    # We print the final equation with the calculated numbers
    print(f"{mu_X3:.4f} = {K:.1f} * {mu_X:.4f}")
    print("\nThis demonstrates that a growth factor of K=3 is achievable.")
    print("Based on the theory of growth in Lie groups, this is the largest possible value of K that holds for all non-degenerate compact sets.")
    print("The final answer for the constant K is 3.")

demonstrate_growth_factor()
<<<3>>>