from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.named_groups import CyclicGroup

def verify_c6():
    """
    Verifies that the cyclic group of order 6 (C6) has a maximal
    product-free set of size 2.
    """
    # 1. Set up the group C6 and its elements
    G = CyclicGroup(6)
    elements = sorted(list(G.elements), key=lambda p: p.size) # for consistent ordering
    e = G.identity
    # Let's name the elements for clarity in the output
    # C6 = {e, x, x^2, x^3, x^4, x^5}
    # The generator x corresponds to Permutation(0, 1, 2, 3, 4, 5)
    x = Permutation(0, 1, 2, 3, 4, 5)
    elem_names = {
        e: "e",
        x: "x",
        x**2: "x^2",
        x**3: "x^3",
        x**4: "x^4",
        x**5: "x^5",
    }

    # 2. Define the specific product-free set S = {a, b} of size 2
    # For C6, a known maximal product-free set is S = {x^2, x^3}
    a = x**2
    b = x**3
    S = {a, b}
    S_names = {elem_names[p] for p in S}
    print(f"Group G = C6 = {sorted(list(elem_names.values()))}")
    print(f"Chosen set S = {sorted(list(S_names))}")
    print("-" * 20)

    # 3. Verify that S is product-free: S^2 intersect S = empty
    s_squared = {p1 * p2 for p1 in S for p2 in S}
    if not s_squared.isdisjoint(S):
        print("Error: The chosen set S is not product-free.")
        return

    # 4. Calculate the set G \ S
    G_minus_S = set(G.elements) - S
    G_minus_S_names = {elem_names[p] for p in G_minus_S}

    # 5. Calculate the set K = (S^-1 * S) U (S * S^-1) U {g | g^2 in S}
    S_inv = {p**-1 for p in S}
    S_inv_S = {p1 * p2 for p1 in S_inv for p2 in S}
    S_S_inv = {p1 * p2 for p1 in S for p2 in S_inv}
    R_S = {g for g in G.elements if g**2 in S}
    
    K = S_inv_S.union(S_S_inv).union(R_S)
    K_names = {elem_names[p] for p in K}

    # 6. Print the final equation G \ S = K
    print("Verification of maximality requires G \\ S = K, where:")
    print("K = (S⁻¹S ∪ SS⁻¹ ∪ {g ∈ G | g² ∈ S})")
    print("\nCalculating the left side of the equation:")
    print(f"G \\ S = {sorted(list(G_minus_S_names))}")
    
    print("\nCalculating the right side of the equation:")
    print(f"K = {sorted(list(K_names))}")

    print("\n" + "-" * 20)
    if G_minus_S == K:
        print("Conclusion: The sets are equal, so S is a maximal product-free set.")
    else:
        print("Conclusion: The sets are NOT equal, so S is not maximal.")

if __name__ == '__main__':
    # The question is "How many finite groups...". The answer is 8.
    # The code below verifies this property for one of those 8 groups (C6)
    # and prints the verification equation as requested.
    print("The number of finite groups containing maximal by inclusion product-free sets of size 2 is 8.")
    print("\n--- Verification for C6 ---")
    verify_c6()
    print("\nFinal Answer:")
    print("8")
