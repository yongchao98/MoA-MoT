import math

def nCr_poly(n, r):
    """
    Computes binomial coefficient C(n,r)
    """
    if r < 0 or r > n:
        return 0
    if r == 0 or r == n:
        return 1
    if r > n // 2:
        r = n - r
    
    res = 1
    for i in range(r):
        res = res * (n - i) // (i + 1)
    return res

def find_closed_form_coefficients():
    """
    This function calculates the coefficients for the closed form of the sum.
    The sum Sn can be expressed as:
    Sn = 4^n * sum_{i=0 to 5} c_i * C(n, i)
    The function calculates these integer coefficients c_i.
    """

    # Coefficients of the polynomial P_5(x) found from the recurrence.
    # P_5(x) = 1 + 464x + 8256x^2 + 18944x^3 + 4096x^4
    coeffs_P5 = [1, 464, 8256, 18944, 4096]

    # The sum is given by Sn = [x^n] P_5(x) / (1-4x)^6
    # This leads to a form Sn = 4^n * sum_{k=0 to 4} a_k * C(n+5-k, 5)
    # where a_k = coeffs_P5[k] / 4^k
    
    a = [coeffs_P5[k] // (4**k) for k in range(5)]
    # a = [1, 116, 516, 296, 16]

    # We convert the basis {C(n+j, 5)} to {C(n, i)}
    # C(n+j, 5) = sum_{i=0 to j} C(j, i) * C(n, 5-i)
    # The final expression is Sn = 4^n * sum_{i=0 to 5} c_i * C(n, i)
    # The coefficients c_i are computed below.
    
    final_coeffs = [0] * 6 # for c_0 to c_5
    
    for i in range(6): # Corresponds to coefficient for C(n, 5-i)
        c_index = 5 - i
        # c_{5-i} = sum_{k=0 to 4} a_k * C(5-k, i)
        current_c = 0
        for k in range(5):
             current_c += a[k] * nCr_poly(5 - k, i)
        final_coeffs[c_index] = current_c
        
    c5, c4, c3, c2, c1, c0 = final_coeffs
    
    print("The closed form for the sum is:")
    print("S_n = 4^n * (c_5*C(n,5) + c_4*C(n,4) + c_3*C(n,3) + c_2*C(n,2) + c_1*C(n,1) + c_0*C(n,0))")
    print("\nWhere C(n,k) is the binomial coefficient, and the integer coefficients are:")
    print(f"c_5 = {c5}")
    print(f"c_4 = {c4}")
    print(f"c_3 = {c3}")
    print(f"c_2 = {c2}")
    print(f"c_1 = {c1}")
    print(f"c_0 = {c0}")

find_closed_form_coefficients()