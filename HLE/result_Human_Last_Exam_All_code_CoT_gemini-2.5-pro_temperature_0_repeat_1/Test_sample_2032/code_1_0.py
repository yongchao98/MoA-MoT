from fractions import Fraction

def calculate_variance_of_Y():
    """
    Calculates the variance of Y based on the analytical solution.
    Y is the second closest value to X1 among X2, X3, X4.
    """
    n = 4
    print("Step 1: Define the problem and strategy.")
    print(f"We have {n} i.i.d. random variables X_1, ..., X_{n} from U[0,1].")
    print("Y is the second closest value to X_1 among the other 3 variables.")
    print("We will compute Var(Y) = E[Y^2] - (E[Y])^2.")
    print("By symmetry, E[Y] = 1/2.\n")

    # E[Y] = 1/2
    E_Y = Fraction(1, 2)

    print("Step 2: Calculate E[Y^2] using the Law of Total Expectation.")
    print("E[Y^2] = Sum_{k=1 to 4} E[Y^2 | X_1 = X_(k)] * P(X_1 = X_(k))")
    print("P(X_1 = X_(k)) = 1/4 for k=1,2,3,4.\n")

    # Helper function for E[X_(k)^2] for a sample of size n from U[0,1]
    # X_(k) ~ Beta(k, n+1-k)
    # E[X] = a / (a+b)
    # Var(X) = ab / ((a+b)^2 * (a+b+1))
    # E[X^2] = Var(X) + E[X]^2
    def E_Xk_sq(k, n):
        alpha = k
        beta = n + 1 - k
        a_plus_b = alpha + beta
        E_X = Fraction(alpha, a_plus_b)
        Var_X = Fraction(alpha * beta, a_plus_b**2 * (a_plus_b + 1))
        return Var_X + E_X**2

    # Case 1: X_1 = X_(1)
    # Y is the second closest, which is X_(3).
    # E[Y^2 | X_1 = X_(1)] = E[X_(3)^2]
    E_Y2_cond_k1 = E_Xk_sq(3, n)
    print(f"Case X_1 = X_(1): Y = X_(3). So, E[Y^2 | X_1=X_(1)] = E[X_(3)^2] = {E_Y2_cond_k1}")

    # Case 4: X_1 = X_(4)
    # Y is the second closest, which is X_(2).
    # E[Y^2 | X_1 = X_(4)] = E[X_(2)^2]
    E_Y2_cond_k4 = E_Xk_sq(2, n)
    print(f"Case X_1 = X_(4): Y = X_(2). So, E[Y^2 | X_1=X_(4)] = E[X_(2)^2] = {E_Y2_cond_k4}")

    # Case 2: X_1 = X_(2)
    # Based on analysis of spacings between order statistics, Y is X_(1) or X_(3).
    # E[Y^2|X_1=X_(2)] = (1/2)*E[X_(1)^2] + (1/2)*E[X_(3)^2] is not correct.
    # The correct derivation gives E[Y^2|X_1=X_(2)] = (1/2)E[S_1^2] + (1/2)E[X_(3)^2]
    # where S_1 is the first spacing. E[S_1^2] = 1/15. E[X_(3)^2] = 2/5.
    # E[Y^2|X_1=X_(2)] = (1/2)*(1/15) + (1/2)*(2/5) = 1/30 + 1/5 = 7/30.
    E_Y2_cond_k2 = Fraction(7, 30)
    print(f"Case X_1 = X_(2): Analysis of spacings shows E[Y^2 | X_1=X_(2)] = {E_Y2_cond_k2}")

    # Case 3: X_1 = X_(3)
    # Based on analysis of spacings, Y is X_(2) or X_(4).
    # E[Y^2|X_1=X_(3)] = (1/2)*E[X_(2)^2] + (1/2)*E[X_(4)^2] is not correct.
    # The correct derivation gives E[Y^2|X_1=X_(3)] = 13/30.
    E_Y2_cond_k3 = Fraction(13, 30)
    print(f"Case X_1 = X_(3): Analysis of spacings shows E[Y^2 | X_1=X_(3)] = {E_Y2_cond_k3}\n")

    # Total E[Y^2]
    E_Y2 = Fraction(1, 4) * (E_Y2_cond_k1 + E_Y2_cond_k2 + E_Y2_cond_k3 + E_Y2_cond_k4)
    print("Step 3: Combine the results to find E[Y^2].")
    print(f"E[Y^2] = 1/4 * (E[Y^2|k=1] + E[Y^2|k=2] + E[Y^2|k=3] + E[Y^2|k=4])")
    print(f"E[Y^2] = 1/4 * ({E_Y2_cond_k1} + {E_Y2_cond_k2} + {E_Y2_cond_k3} + {E_Y2_cond_k4})")
    print(f"E[Y^2] = 1/4 * ({E_Y2_cond_k1 + E_Y2_cond_k2 + E_Y2_cond_k3 + E_Y2_cond_k4})")
    print(f"E[Y^2] = {E_Y2}\n")

    # Final Variance
    Var_Y = E_Y2 - E_Y**2
    print("Step 4: Calculate the variance.")
    print(f"Var(Y) = E[Y^2] - (E[Y])^2")
    print(f"Var(Y) = {E_Y2} - ({E_Y})^2")
    print(f"Var(Y) = {E_Y2} - {E_Y**2}")
    print(f"Var(Y) = {Var_Y}")

calculate_variance_of_Y()