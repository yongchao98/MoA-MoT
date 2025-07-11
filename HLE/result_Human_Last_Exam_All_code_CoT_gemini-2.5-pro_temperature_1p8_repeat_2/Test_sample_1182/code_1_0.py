import math

def solve_curve_reduction():
    """
    This script finds the stable reduction of a hyperelliptic curve
    and determines the number of double points.
    
    The curve is y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5
    We work over the prime p = 2.
    """
    
    # Coefficients of f(x) for powers x^1, x^2, x^3, x^4, x^5
    c_orig = [8, 1, 4, 4, 8]
    p = 2
    
    # Through analysis, we find the best transformation uses e=3.
    # The transformation is x = p^e * X, y = p^k * Y
    e = 3
    
    # We choose k such that the coefficient of Y^2 becomes 1.
    # The condition is that 2k equals the minimum valuation of the terms on the RHS.
    # min_valuation = min(v_p(c_i) + i*e for i in range(1, 6))
    # v_p(8) = 3, v_p(1) = 0, v_p(4) = 2, v_p(4) = 2, v_p(8) = 3
    v = [3, 0, 2, 2, 3]
    valuations = [v[i-1] + i*e for i in range(1, 6)]
    min_valuation = min(valuations) # This is 6 for e=3
    
    # We need 2k = min_valuation
    k = min_valuation // 2

    print(f"Starting with the curve: y^2 = {c_orig[0]}*x + {c_orig[1]}*x^2 + {c_orig[2]}*x^3 + {c_orig[3]}*x^4 + {c_orig[4]}*x^5")
    print(f"To find the stable reduction, we perform a change of variables: x = {p**e}*X and y = {p**k}*Y.")
    print(f"This is x = {int(p**e)}*X and y = {int(p**k)}*Y.")
    
    # The transformed equation is Y^2 = sum(c'_i * X^i), where c'_i = c_i * p^(i*e - 2k)
    new_c = []
    for i in range(1, 6):
        coeff = c_orig[i-1] * (p**(i*e - 2*k))
        new_c.append(int(coeff))

    print("\nThe transformed equation is Y^2 = c'_1*X + c'_2*X^2 + c'_3*X^3 + c'_4*X^4 + c'_5*X^5.")
    print("The new integer coefficients are:")
    for i in range(5):
        print(f"c'_{i+1} = {new_c[i]}")
    
    print("\nThe full transformed equation is:")
    final_equation_str = f"Y^2 = {new_c[0]}*X + {new_c[1]}*X^2 + {new_c[2]}*X^3 + {new_c[3]}*X^4 + {new_c[4]}*X^5"
    print(final_equation_str)

    print("\nReducing this equation modulo 2, we get:")
    reduced_parts = []
    for i, coeff in enumerate(new_c):
        if coeff % 2 != 0:
            power = i + 1
            if power == 1:
                reduced_parts.append("X")
            else:
                reduced_parts.append(f"X^{power}")
    reduced_equation = "Y^2 = " + " + ".join(reduced_parts)
    print(reduced_equation)
    
    print("\nThe curve Y^2 = X + X^2 is smooth, meaning it has no singular points.")
    print("Since double points are a type of singular point, this curve has no double points.")
    
solve_curve_reduction()