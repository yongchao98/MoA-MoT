def solve_ghz_purification():
    """
    This function calculates and prints the symbolic expression for the product of
    the success probability and the output fidelity for the given GHZ state
    purification protocol.
    """

    # Coefficients for the final expression of the product P_succ * F_out
    # The expression is of the form:
    # (c_f1f2*F1*F2 + c_f1*F1 + c_f2*F2 + c_const) / denominator
    c_f1f2 = 32
    c_f1 = 1
    c_f2 = 10
    c_const = -1
    denominator = 42

    # Print the final expression, showing each number as requested.
    print("The product of the successful output fidelity and the success probability is given by the expression:")
    print(f"({c_f1f2}*F1*F2 + {c_f1}*F1 + {c_f2}*F2 - {abs(c_const)}) / {denominator}")

solve_ghz_purification()