def print_closed_form():
    """
    This function prints the derived closed form for the sum:
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)
    
    The closed form is of the structure: 4^n * P(n),
    where P(n) is a polynomial in n of degree 5.
    We express P(n) in the basis of binomial coefficients C(n, k).
    """
    
    # The coefficients of the polynomial in the binomial basis C(n, k)
    # for k from 5 down to 0 have been calculated.
    coeffs = [973, 2681, 2578, 990, 121, 1]
    
    # Construct and print the formula string.
    # The notation C(n, k) is used for the binomial coefficient "n choose k".
    formula = (
        f"4**n * ({coeffs[0]}*C(n, 5) + {coeffs[1]}*C(n, 4) + "
        f"{coeffs[2]}*C(n, 3) + {coeffs[3]}*C(n, 2) + "
        f"{coeffs[4]}*C(n, 1) + {coeffs[5]}*C(n, 0))"
    )
    
    print("The closed form for the sum is:")
    print(formula)

print_closed_form()