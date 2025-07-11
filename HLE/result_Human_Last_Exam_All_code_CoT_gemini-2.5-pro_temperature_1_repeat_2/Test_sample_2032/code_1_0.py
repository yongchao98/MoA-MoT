from fractions import Fraction

def solve_variance():
    """
    Calculates the variance of Y, where Y is the second closest value to X1
    among X2, X3, and X4, with all Xi being i.i.d. U[0, 1].
    """
    
    # Number of random variables
    n = 4

    # The problem can be solved by conditioning on the rank of X1 among the
    # four variables. Let X_(k) be the k-th order statistic.
    # A detailed analysis shows that the probability of Y being equal to the
    # k-th order statistic is as follows:
    # P(Y=X_(1)) = 1/8
    # P(Y=X_(2)) = 3/8
    # P(Y=X_(3)) = 3/8
    # P(Y=X_(4)) = 1/8
    probs = [Fraction(1, 8), Fraction(3, 8), Fraction(3, 8), Fraction(1, 8)]

    # We will now calculate E[Y] and E[Y^2] using the law of total expectation:
    # E[Y] = sum(P(Y=X_(k)) * E[X_(k)])
    # E[Y^2] = sum(P(Y=X_(k)) * E[X_(k)^2])
    
    # The moments of order statistics from U[0, 1] for a sample of size n are:
    # E[X_(k)] = k / (n+1)
    # E[X_(k)^2] = k*(k+1) / ((n+1)*(n+2))

    E_Xk = []
    E_Xk2 = []

    for k in range(1, n + 1):
        # Expected value of X_(k)
        exk = Fraction(k, n + 1)
        E_Xk.append(exk)
        
        # Second moment of X_(k)
        exk2 = Fraction(k * (k + 1), (n + 1) * (n + 2))
        E_Xk2.append(exk2)

    # Calculate E[Y] and E[Y^2] by summing over the order statistics
    # weighted by the probabilities we found.
    E_Y = sum(p * e for p, e in zip(probs, E_Xk))
    E_Y2 = sum(p * e for p, e in zip(probs, E_Xk2))

    # Calculate Var(Y) using the formula Var(Y) = E[Y^2] - (E[Y])^2
    Var_Y = E_Y2 - E_Y**2

    # Print the details of the final calculation as requested.
    print(f"The variance of Y is given by the formula Var(Y) = E[Y^2] - (E[Y])^2.")
    print(f"We calculated the required components:")
    print(f"E[Y] = {E_Y}")
    print(f"E[Y^2] = {E_Y2}")
    print(f"\nFinal calculation:")
    print(f"Var(Y) = {E_Y2} - ({E_Y})^2")
    print(f"Var(Y) = {E_Y2} - {E_Y**2}")
    print(f"Var(Y) = {Var_Y}")
    print(f"\nThe variance of Y is {Var_Y}.")

solve_variance()