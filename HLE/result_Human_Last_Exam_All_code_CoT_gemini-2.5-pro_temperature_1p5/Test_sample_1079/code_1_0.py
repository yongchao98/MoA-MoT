import math

def solve_h3_regular_elements():
    """
    Calculates the number of elements in the H3 reflection group that have a regular
    eigenvector with a corresponding eigenvalue of order 10.
    """
    # Step 1: Define properties of the H3 reflection group.
    group_name = "H3"
    group_order = 120
    coxeter_number = 10

    print(f"Step 1: The properties of the reflection group {group_name} are:")
    print(f"  - Group Order |{group_name}| = {group_order}")
    print(f"  - Coxeter Number h = {coxeter_number}\n")

    # Step 2: Relate the problem to regular elements.
    print("Step 2: We need to count elements 'w' with a regular eigenvalue of order 10.")
    print("This is equivalent to counting elements that are regular and have order 10.")
    print("These are elements conjugate to c^k, where 'c' is a Coxeter element and gcd(k, h) = 1.\n")

    # Step 3: Identify the distinct conjugacy classes.
    coprime_k = [k for k in range(1, coxeter_number) if math.gcd(k, coxeter_number) == 1]
    num_classes = len(coprime_k) // 2

    print(f"Step 3: The values of k coprime to h={coxeter_number} are {coprime_k}.")
    print(f"For H3, these fall into {num_classes} distinct conjugacy classes (corresponding to {1, 9} and {3, 7}).\n")

    # Step 4: Calculate the size of each conjugacy class.
    class_size = group_order // coxeter_number

    print("Step 4: The size of each such class is |H3| / h.")
    print(f"  - Size = {group_order} / {coxeter_number} = {class_size}\n")
    
    # Step 5: Sum the sizes of the disjoint classes to get the total.
    class_1_size = class_size
    class_2_size = class_size
    total_elements = class_1_size + class_2_size

    print("Step 5: The total number of elements is the sum of the sizes of these classes.")
    print(f"The final calculation is:")
    print(f"Total = {class_1_size} + {class_2_size} = {total_elements}")

solve_h3_regular_elements()