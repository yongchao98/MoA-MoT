def solve_tori_number():
    """
    Calculates the exact number of Fq-rational maximal tori for a reductive group of type E8.
    """
    # The rank of the E8 group
    n = 8

    # The order of the Weyl group of E8 is 2^14 * 3^5 * 5^2 * 7
    W_order = (2**14) * (3**5) * (5**2) * 7

    # For E8, the dimension of the invariant subspace of the k-th exterior power
    # of the reflection representation, a_k, is 1 for k=0 and 0 for k=1,...,8.
    a = [0] * (n + 1)
    a[0] = 1

    # The total number of tori is given by the formula:
    # N = |W| * sum_{k=0 to n} (-1)^k * a_k * q^(n-k)
    # Since only a[0] is non-zero, this simplifies to |W| * a[0] * q^n.

    # We are asked to output each number in the final equation.
    # The final equation is N = W_order * (a[0]*q^8 - a[1]*q^7 + ... + a[8])
    # The term in the parenthesis is (1*q^8 - 0*q^7 + ... + 0) = q^8
    
    print(f"The number of Fq-rational maximal tori is given by the formula: |W(E8)| * (a_0*q^8 - a_1*q^7 + ... + a_8)")
    print(f"The order of the Weyl group W(E8) is {W_order}.")
    print("The coefficients a_k = dim((wedge^k V)^W) are:")
    print(f"a_0 = {a[0]}")
    for k in range(1, n + 1):
        print(f"a_{k} = {a[k]}")
    
    print("\nSubstituting these values, the expression simplifies:")
    equation = f"{W_order} * ({a[0]}*q**8)"
    print(f"N = {equation}")

solve_tori_number()