import math

def solve():
    """
    Solves the group theory problem by calculating the sum of orders of outer automorphism groups.
    """
    # Step 1: Analyze the group G.
    # The group G is given by the presentation:
    # G = <a, b, c, d | aba^-1 = a^2, bcb^-1 = c^2, cdc^-1 = d^2, dad^-1 = a^2>
    # It is a known, though non-trivial, fact of combinatorial group theory that this presentation
    # defines the trivial group, G = {1}. We proceed based on this fact.

    # Step 2: Characterize the central extensions.
    # The set E of isomorphism classes of central extensions of G by an abelian group C
    # is in one-to-one correspondence with Hom(H_2(G, Z), C), where H_2(G, Z) is the Schur multiplier of G.
    # The group C is the cyclic group of order 31.
    order_of_C = 31

    # Step 3: Count the number of extensions.
    # Since G is the trivial group {1}, its Schur multiplier H_2(G, Z) is also trivial.
    # The set of extensions E corresponds to Hom({1}, C_31).
    # There is only one homomorphism from the trivial group to any group, so there is only one
    # central extension E up to isomorphism.

    # Step 4: Identify the unique extension E.
    # This unique extension is the direct product E = G x C.
    # Since G = {1} and C = C_31, we have E is isomorphic to C_31.

    # Step 5: Compute the order of the outer automorphism group of E.
    # We need to find the order of Out(E) = Out(C_31).
    # The outer automorphism group is Out(H) = Aut(H) / Inn(H).
    # For an abelian group like C_31, the inner automorphism group Inn(C_31) is trivial.
    # Therefore, Out(C_31) is isomorphic to Aut(C_31).

    # Step 6: Compute the order of the automorphism group of C_31.
    # The order of the automorphism group of the cyclic group C_n is given by Euler's totient function, phi(n).
    # We need to compute |Aut(C_31)| = phi(31).
    # Since 31 is a prime number, phi(31) = 31 - 1.
    order_of_outer_automorphism_group = order_of_C - 1

    # Step 7: Compute the final sum.
    # The collection E contains only one element, E.
    # The sum is therefore just the order of the outer automorphism group of E.
    final_sum = order_of_outer_automorphism_group

    print("The problem asks for the sum of orders of outer automorphism groups of all central extensions of G by C.")
    print("1. The group G is the trivial group {1}.")
    print("2. There is only one central extension E, which is isomorphic to G x C = C_31.")
    print("3. The order of the outer automorphism group of E is o(E) = |Out(C_31)|.")
    print("4. Since C_31 is abelian, |Out(C_31)| = |Aut(C_31)|.")
    print("5. The order of Aut(C_n) is given by Euler's totient function phi(n).")
    print("\nFinal Calculation:")
    print(f"Sum = |Aut(C_{order_of_C})| = phi({order_of_C}) = {order_of_C} - 1 = {final_sum}")

solve()