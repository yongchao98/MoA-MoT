import math

def solve_group_theory_problem():
    """
    Solves the group theory problem by analyzing the groups and computing the sum.
    """

    # The group C is the cyclic group of order 31.
    order_of_C = 31

    # Plan Step 1 & 2: Analyze group G and the extensions.
    # The group G is given by the presentation:
    # <a, b, c, d | aba^{-1} = a^2, bcb^{-1} = c^2, cdc^{-1} = d^2, dad^{-1} = a^2>
    # From the first relation, aba^{-1} = a^2, we can deduce b = a^2.
    # Substituting b = a^2 into the second relation gives a^2ca^{-2} = c^2.
    # The relations for G (on generators a,c,d) become:
    # 1. dad^{-1} = a^2
    # 2. cdc^{-1} = d^2
    # 3. a^2ca^{-2} = c^2
    # From these relations, one can show that all generators a,c,d lie in the
    # commutator subgroup [G,G], making G a perfect group. Such groups with these
    # types of "squaring" relations are known to be trivial.
    # Thus, we conclude G = {1}.

    # Plan Step 3 & 4: Characterize the unique extension E.
    # The set of central extensions E is classified by H^2(G, C).
    # Since G = {1}, H^2({1}, C) is the trivial group {0}.
    # This means there is only one central extension, E, up to isomorphism.
    # This extension is the direct product E = G x C = {1} x C, which is isomorphic to C.

    # Plan Step 5, 6, 7: Compute the sum.
    # The sum is over a single element E, so we just need to compute o(E) = o(C).
    # o(C) = |Out(C)| = |Aut(C)/Inn(C)|.
    # C is Z_31, an abelian group, so its inner automorphism group Inn(C) is trivial (|Inn(C)| = 1).
    # Therefore, o(C) = |Aut(C)|.
    # The order of the automorphism group of the cyclic group Z_n is given by Euler's totient function, phi(n).
    
    # We need to compute phi(31).
    p = order_of_C
    
    # Since 31 is a prime number, phi(31) = 31 - 1.
    order_of_aut_c = p - 1
    
    # The sum is just this value.
    the_sum = order_of_aut_c

    print(f"The group C is the cyclic group of order {p}.")
    print("Based on analysis of the presentation, the group G is the trivial group {1}.")
    print("Therefore, there is only one central extension E, which is isomorphic to C.")
    print("The problem reduces to finding the order of the outer automorphism group of C, o(C).")
    print(f"o(C) = |Aut(Z_{p})| because C is abelian.")
    print(f"The order of the automorphism group of Z_{p} is phi({p}).")
    print("The final equation is:")
    print(f"Sum = phi({p}) = {p} - 1 = {the_sum}")

solve_group_theory_problem()