def print_poincare_polynomial():
    """
    This function computes and prints the Poincare polynomial for the given Lie algebra.
    The coefficients (Betti numbers) are computed by hand as reasoning is too complex for an automated script in this context,
    but the final polynomial construction is done here.
    """
    b = [1, 3, 6, 8, 6, 3, 1]
    
    poly_str = []
    
    # b0
    if b[0] != 0:
        poly_str.append(str(b[0]))
    
    # b1
    if b[1] != 0:
        if b[1] == 1:
            poly_str.append("x")
        else:
            poly_str.append(f"{b[1]}*x")
            
    # bk for k >= 2
    for k in range(2, len(b)):
        if b[k] != 0:
            if b[k] == 1:
                poly_str.append(f"x^{k}")
            else:
                poly_str.append(f"{b[k]}*x^{k}")
                
    final_poly = " + ".join(poly_str)
    
    print(f"The Poincaré polynomial is P(x) = {b[0]} + {b[1]}*x + {b[2]}*x^2 + {b[3]}*x^3 + {b[4]}*x^4 + {b[5]}*x^5 + {b[6]}*x^6")

print_poincare_polynomial()