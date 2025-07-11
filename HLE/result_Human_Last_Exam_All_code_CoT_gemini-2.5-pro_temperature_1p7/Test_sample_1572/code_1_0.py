def solve_p_n_formula():
    """
    This function derives and prints the formula for P(n).
    """
    # The coefficients c_2, c_4 for the expansion of ln(Q(n)/T(n)) are:
    # ln(Q(n)/T(n)) = c_2/n^2 + c_4/n^4 + O(n^-6)
    c2_num = 1
    c2_den = 720
    c4_num = -1
    c4_den = 5040
    
    # To achieve a relative error of O(n^-6), ln(P(n)) must match these terms.
    # ln(P(n)) = c_2/n^2 + c_4/n^4
    # Therefore, P(n) is the exponential of this expression.
    
    print("The formula for P(n) is:")
    # Using python to format the string with numbers to clearly show the result.
    formula = "P(n) = exp(1/({0} * n^2) - 1/({1} * n^4))".format(c2_den, c4_den)
    print(formula)

solve_p_n_formula()