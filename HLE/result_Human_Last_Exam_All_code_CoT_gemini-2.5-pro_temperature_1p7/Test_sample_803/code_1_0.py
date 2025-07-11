import sympy

def investigate_filled_groups(q, m):
    """
    Investigates whether nonabelian groups of order 2*q**m can be filled,
    based on a key theorem from group theory.
    
    The code walks through a mathematical proof to arrive at the answer.
    """
    
    # Input validation
    try:
        if not (isinstance(q, int) and isinstance(m, int) and m >= 1 and q > 2 and sympy.isprime(q)):
             raise ValueError("InputError")
    except (TypeError, ValueError):
        print("Error: Please provide an odd prime number for q (e.g., 3, 5, 7) and a natural number (>=1) for m.")
        return

    n = q**m
    order = 2 * n

    print(f"Investigating filled groups of order 2 * q^m = 2 * {q}^{m} = {order}")
    print("="*60)

    # Step 1: State the relevant theorem
    print("A key theorem states a finite group G is 'filled' if and only if for EVERY maximal subgroup M,")
    print("the quotient group G/core_G(M) is NOT a cyclic group of a prime order.")
    print("This means G is NOT filled if we can find just ONE maximal subgroup M")
    print("for which G/core_G(M) IS a cyclic group of a prime order.\n")

    # Step 2: Analyze the structure of any group G of the given order
    print(f"Any group G of order {order} must have a normal subgroup of index 2.")
    print("\n--- Mathematical Argument ---")
    print(f"1. By Sylow's Theorems, the number of Sylow q-subgroups (each of order q^m = {n})")
    print(f"   must divide 2 and also be congruent to 1 (mod {q}). The only possibility is 1.")
    print(f"2. Therefore, G has a unique Sylow q-subgroup, let's call it Q.")
    print(f"3. A unique Sylow subgroup is always a normal subgroup in G.")
    print(f"4. The index of this normal subgroup Q in G is |G|/|Q| = {order} / {n} = 2.\n")
    
    # Step 3: Apply the theorem to this structure
    print("--- Applying the Theorem ---")
    print(f"1. A subgroup of index 2 is always a maximal subgroup. So, Q is a maximal subgroup of G.")
    print(f"2. The 'core' of a subgroup M, core_G(M), is the largest normal subgroup of G contained within M.")
    print(f"3. Since Q is itself a normal subgroup, core_G(Q) = Q.")
    print(f"4. We must check the quotient group for this maximal subgroup: G/core_G(Q), which simplifies to G/Q.")
    print(f"5. The order of this quotient group is |G|/|Q| = {order} / {n} = 2.")
    print(f"6. Any group of prime order is cyclic. So, G/Q is a cyclic group of prime order 2 (isomorphic to C_2).\n")

    # Step 4: Conclusion
    print("--- Conclusion ---")
    print("We have found a maximal subgroup (Q) for which G/core_G(Q) IS a cyclic group of prime order.")
    print("According to the theorem, this is sufficient to prove that the group G is NOT filled.")
    print("This logic is valid for ANY group of order 2*q^m, including all nonabelian ones.")
    
    # Final 'equation' breakdown as requested.
    print("\nFinal analysis based on the disproving condition:")
    final_equation_index = 2
    final_equation_order = 2
    print(f"Index of maximal subgroup Q is [G:Q] = {final_equation_index}")
    print(f"Order of quotient G/core_G(Q) is |G/Q| = {final_equation_order}")
    print(f"Is the order {final_equation_order} a prime number? {sympy.isprime(final_equation_order)}")
    print("The condition for a group to be 'filled' is therefore violated.")

# Example call. You can change q and m to any odd prime and natural number.
# Let's try q=5, m=1. This corresponds to groups of order 10, e.g., the Dihedral group D_10.
investigate_filled_groups(q=5, m=1)