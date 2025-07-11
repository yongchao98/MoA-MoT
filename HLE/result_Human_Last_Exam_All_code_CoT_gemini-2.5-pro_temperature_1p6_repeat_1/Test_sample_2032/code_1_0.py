from fractions import Fraction

def beta_moment2(alpha, beta):
    """Calculates the second moment E[X^2] for a Beta(alpha, beta) variable."""
    return (alpha * (alpha + 1)) / ((alpha + beta) * (alpha + beta + 1))

def main():
    """
    Calculates the variance of Y based on the analytical solution.
    """
    n = 4
    print("Step 1: Calculate E[Y]")
    # By symmetry E[Y] = 1/2
    E_Y = Fraction(1, 2)
    print(f"By symmetry, E[Y] = {E_Y}\n")

    print("Step 2: Calculate the second moments of the order statistics X_(k) for n=4.")
    E_X_k_sq = {}
    for k in range(1, n + 1):
        alpha = k
        beta = n - k + 1
        moment = beta_moment2(alpha, beta)
        E_X_k_sq[k] = Fraction(moment).limit_denominator()
        print(f"E[X_({k})^2] = {E_X_k_sq[k].numerator}/{E_X_k_sq[k].denominator}")
    
    print("\nStep 3: Calculate the conditional expectations E[Y^2 | X_1 = X_(k)].")
    E_Y_sq_cond = {}
    
    # Case X_1 = X_(1) => Y = X_(3)
    E_Y_sq_cond[1] = E_X_k_sq[3]
    print(f"E[Y^2 | X_1=X_(1)] = E[X_(3)^2] = {E_Y_sq_cond[1].numerator}/{E_Y_sq_cond[1].denominator}")

    # Case X_1 = X_(4) => Y = X_(2)
    E_Y_sq_cond[4] = E_X_k_sq[2]
    print(f"E[Y^2 | X_1=X_(4)] = E[X_(2)^2] = {E_Y_sq_cond[4].numerator}/{E_Y_sq_cond[4].denominator}")
    
    # Case X_1 = X_(2)
    E_Y_sq_cond[2] = (E_X_k_sq[1] + E_X_k_sq[3]) / 2
    val1, val3 = E_X_k_sq[1], E_X_k_sq[3]
    print(f"E[Y^2 | X_1=X_(2)] = 1/2 * (E[X_(1)^2] + E[X_(3)^2]) = 1/2 * ({val1.numerator}/{val1.denominator} + {val3.numerator}/{val3.denominator}) = {E_Y_sq_cond[2].numerator}/{E_Y_sq_cond[2].denominator}")

    # Case X_1 = X_(3)
    E_Y_sq_cond[3] = (E_X_k_sq[2] + E_X_k_sq[4]) / 2
    val2, val4 = E_X_k_sq[2], E_X_k_sq[4]
    print(f"E[Y^2 | X_1=X_(3)] = 1/2 * (E[X_(2)^2] + E[X_(4)^2]) = 1/2 * ({val2.numerator}/{val2.denominator} + {val4.numerator}/{val4.denominator}) = {E_Y_sq_cond[3].numerator}/{E_Y_sq_cond[3].denominator}")

    print("\nStep 4: Calculate E[Y^2] using the law of total expectation.")
    E_Y_sq = sum(E_Y_sq_cond.values()) / n
    terms_str = [f"{v.numerator}/{v.denominator}" for v in E_Y_sq_cond.values()]
    print(f"E[Y^2] = 1/4 * ({' + '.join(terms_str)}) = {E_Y_sq.numerator}/{E_Y_sq.denominator}\n")

    print("Step 5: Calculate Var(Y).")
    Var_Y = E_Y_sq - E_Y**2
    e_y_sq_n, e_y_sq_d = E_Y_sq.numerator, E_Y_sq.denominator
    e_y_n, e_y_d = E_Y.numerator, E_Y.denominator
    print(f"Var(Y) = E[Y^2] - (E[Y])^2 = {e_y_sq_n}/{e_y_sq_d} - ({e_y_n}/{e_y_d})^2 = {Var_Y.numerator}/{Var_Y.denominator}")
    
    final_answer = float(Var_Y)
    print(f"\nThe final variance is {final_answer:.6f}")


if __name__ == '__main__':
    main()
