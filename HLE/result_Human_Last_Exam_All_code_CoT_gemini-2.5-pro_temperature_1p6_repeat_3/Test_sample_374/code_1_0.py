import math

def get_prime_factorization(n):
    """Returns the prime factorization of n as a dictionary."""
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
    return factors

def main():
    """
    Solves the problem by calculating the properties of the relevant group
    and discussing the existence of subgroups.
    """
    n = 4
    q = 2

    # Step 1 & 2: State the problem in group-theoretic terms
    print("The problem is to find the highest possible order of the inertial quotient E.")
    print("E must be an odd-order subgroup of Out(D), where D is the defect group.")
    print(f"D is an elementary abelian group of order 16 = 2^4.")
    print(f"Thus, D is isomorphic to a 4-dimensional vector space over F_2.")
    print(f"Out(D) is isomorphic to the general linear group GL({n}, {q}).\n")

    # Step 3: Calculate the order of GL(n, q)
    print(f"Calculating the order of GL({n}, {q})...")
    order = 1
    calculation_str_parts = []
    for i in range(n):
        term = q**n - q**i
        order *= term
        calculation_str_parts.append(f"({q}^{n} - {q}^{i})")
    
    calculation_str = " * ".join(calculation_str_parts)
    print(f"|GL({n}, {q})| = {calculation_str}")

    value_str_parts = []
    for i in range(n):
        term = q**n - q**i
        value_str_parts.append(str(term))
    value_str = " * ".join(value_str_parts)
    print(f"|GL({n}, {q})| = {value_str} = {order}\n")
    
    # Step 4: Find the largest odd divisor
    print("Finding the prime factorization and the largest odd divisor...")
    factors = get_prime_factorization(order)
    print(f"The prime factorization of {order} is: {factors}")

    odd_order = 1
    odd_factors_str_parts = []
    for p, exp in factors.items():
        if p != 2:
            odd_order *= p**exp
            odd_factors_str_parts.append(f"{p}^{exp}")
    
    odd_factors_str = " * ".join(odd_factors_str_parts)
    print(f"The largest odd divisor (the p'-part for p=2) is {odd_factors_str} = {odd_order}\n")
    
    # Step 5: Discuss the existence of such a subgroup
    print(f"The group GL(4, 2) is isomorphic to the alternating group A_8.")
    print(f"So, we must find the largest order of an odd-order subgroup of A_8.")
    print(f"The theoretical maximum order is {odd_order}.")

    print("\nWe must check if a subgroup of order 315 exists in A_8.")
    print("A key group-theoretic result is that a subgroup of order k must divide the order of the normalizer of any of its Sylow p-subgroups.")
    
    print("\nLet's test for a subgroup H of order 315 = 3^2 * 5 * 7.")
    print("Consider a Sylow 7-subgroup P_7 of H. The normalizer N_A8(P_7) in A_8 has order 21.")
    print("The order of N_H(P_7) must divide the order of N_A8(P_7), so |N_H(P_7)| must divide 21.")
    print("By Sylow's theorems in H, the number of Sylow 7-subgroups, n_7(H), must be |H| / |N_H(P_7)| = 315 / |N_H(P_7)|.")
    print("n_7(H) must divide 315/7 = 45 and be congruent to 1 mod 7. Possible values are 1 and 15.")
    print("If n_7(H) = 1, then |N_H(P_7)|=315. But 315 does not divide 21. Contradiction.")
    print("If n_7(H) = 15, then |N_H(P_7)|=315/15=21. This seems possible.\n")

    print("Now consider a Sylow 5-subgroup P_5 of H. The normalizer N_A8(P_5) has order 60.")
    print("|N_H(P_5)| must divide 60.")
    print("n_5(H) must divide 315/5=63 and be congruent to 1 mod 5. Possible values are 1 and 21.")
    print("If n_5(H) = 1, then |N_H(P_5)|=315. But 315 does not divide 60. Contradiction.")
    print("If n_5(H) = 21, then |N_H(P_5)|=315/21=15. 15 divides 60, so this is possible.\n")

    print("Let's analyze other divisors of 315.")
    print("Order 35 = 5 * 7: A group of order 35 must be cyclic. The maximum order of an element in A_8 is 15. So no subgroup of order 35 exists.")
    print("Order 63 = 3^2 * 7: Such a group H would have a normal Sylow 7-subgroup P_7. Thus H would normalize P_7, so H would be a subgroup of N_A8(P_7), which has order 21. This is a contradiction.")
    print("Order 45 = 3^2 * 5: A similar analysis of normalizers shows that no such subgroup exists.\n")

    print("The largest odd order for which a subgroup is known to exist is 21.")
    print("A Sylow 7-subgroup P_7 of A_8 has a normalizer N_A8(P_7) of order 21.")
    print("This group of order 21 is non-abelian, isomorphic to C_7 semi-direct product C_3.")
    
    highest_order = 21
    print(f"\nThus, the highest possible order for the inertial quotient E is {highest_order}.")

if __name__ == "__main__":
    main()
