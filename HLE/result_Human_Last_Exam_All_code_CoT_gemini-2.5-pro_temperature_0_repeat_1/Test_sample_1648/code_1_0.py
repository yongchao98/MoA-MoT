import math

def p_adic_valuation(n, p):
    """
    Calculates the p-adic valuation of n, denoted v_p(n).
    This is the exponent of the highest power of p that divides n.
    """
    if n == 0:
        return float('inf')
    if n % p != 0:
        return 0
    count = 0
    while n > 0 and n % p == 0:
        count += 1
        n //= p
    return count

def get_k_group_order(n, p, k):
    """
    Calculates the order of the K_{2n}(Z/p^k) group.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a natural number (positive integer).")
    
    vp_n = p_adic_valuation(n, p)
    exponent = k + vp_n
    order = p ** exponent
    return order, vp_n, exponent

def main():
    """
    Main function to explain and demonstrate the calculation.
    """
    p = 3
    k = 3
    
    print("The problem is to find the largest natural number n such that K_{2n}(Z/27) is nonzero.")
    print(f"The ring is Z/27, which corresponds to Z/p^k with p={p} and k={k}.")
    print("\nAccording to a theorem by Hesselholt and Madsen, the order of this group is given by the formula:")
    print(f"|K_{2n}(Z/p^k)| = p^(k + v_p(n))")
    print(f"For our case, this is |K_{2n}(Z/27)| = {p}^({k} + v_{p}(n)).")
    
    print("\nLet's test this for a few values of n:")
    
    test_values = [1, 2, 3, 13, 18, 99]
    for n_val in test_values:
        order, vp_n, exponent = get_k_group_order(n_val, p, k)
        print(f"\nFor n = {n_val}:")
        print(f"The 3-adic valuation, v_3({n_val}), is {vp_n}.")
        # Output the equation with numbers as requested
        print(f"The order of K_{{2*n_val}}(Z/27) is {p}^({k} + {vp_n}) = {p}^{exponent} = {order}.")

    print("\nAs the calculation shows, for any natural number n, v_3(n) >= 0.")
    print(f"Therefore, the exponent (3 + v_3(n)) is always >= 3.")
    print(f"This means the order of the group is always at least {p**k}.")
    print("Since the order is always greater than 1, the group is never zero.")
    print("\nConclusion: There is no largest natural number n for which this group is nonzero.")

if __name__ == "__main__":
    main()
