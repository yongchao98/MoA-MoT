def solve_and_print_normal_cone():
    """
    This function provides an explicit representation of the normal cone T_F^°(x^*).

    The normal cone T_F^°(x^*) is the set of vectors s = (s1, s2, s3)
    in R^3 satisfying a specific linear inequality. Based on our analysis,
    this inequality is s3 >= 0.

    To adhere to the output format, we express this as:
    a1*s1 + a2*s2 + a3*s3 >= b
    """

    # Coefficients for the inequality a1*s1 + a2*s2 + a3*s3 >= b
    a1 = 0
    a2 = 0
    a3 = 1
    b = 0

    print("An explicit representation of the normal cone T_F^°(x^*) is the set of all vectors s = (s1, s2, s3) in R^3 that satisfy the inequality:")
    print(f"({a1}) * s1 + ({a2}) * s2 + ({a3}) * s3 >= {b}")

solve_and_print_normal_cone()