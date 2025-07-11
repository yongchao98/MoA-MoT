import math

def phi(n):
    """
    Computes Euler's totient function phi(n).
    """
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def solve_e8_torsion_problem():
    """
    Calculates the number of specified torsion elements in A(E8)/Z.
    """
    # The order of the torsion elements we are interested in.
    d = 10

    # The degrees of the fundamental invariants of the Weyl group W(E_8).
    degrees_E8 = [2, 8, 12, 14, 18, 20, 24, 30]

    # The Coxeter number h is the largest degree, h=30.
    # The minimal word length of the elements is h*(1-1/d) = 30*(9/10) = 27.

    # The number of elements is equal to the number of regular elements of order d
    # in the corresponding Coxeter group W(E_8).

    # Step 1: Calculate the order of W(E_8), which is the product of its degrees.
    order_W_E8 = 1
    for deg in degrees_E8:
        order_W_E8 *= deg

    # Step 2: Calculate phi(d), the Euler's totient function of d.
    phi_d = phi(d)

    # Step 3: Calculate the total number of regular elements of order d.
    # The formula is N = phi(d) * |W| / d.
    num_elements = (phi_d * order_W_E8) // d

    # Print the explanation and the final equation.
    print("The number of requested elements corresponds to the number of regular elements of order 10 in the Weyl group W(E_8).")
    print("This is calculated using the formula: N = phi(d) * |W| / d\n")

    print(f"1. The order of the elements is d = {d}.")
    
    print(f"2. The order of the group W(E_8), |W|, is the product of its degrees {degrees_E8}.")
    order_W_E8_str = f"|W(E_8)| = {' * '.join(map(str, degrees_E8))} = {order_W_E8:,}"
    print(order_W_E8_str)

    print(f"\n3. Euler's totient function for d={d} is phi({d}) = {phi_d}.")

    print("\nFinal Calculation:")
    print(f"Number of elements = phi({d}) * |W(E_8)| / {d}")
    print(f"                     = {phi_d} * {order_W_E8} / {d}")
    print(f"                     = {num_elements}")

solve_e8_torsion_problem()