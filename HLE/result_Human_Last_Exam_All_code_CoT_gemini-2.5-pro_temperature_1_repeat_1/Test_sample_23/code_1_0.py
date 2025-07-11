def solve_admissible_integers():
    """
    Calculates the number of non-admissible integers based on matrix dimensions a and b.

    An integer k is "admissible" if there exist complex a x b matrices A_1,...,A_{ab}
    satisfying:
    1. Each A_i is nonzero
    2. tr(A_i^\dagger A_j) = 0 whenever i != j
    3. Exactly k of the matrices A_i have rank 1.

    This function asks for positive integers a and b and prints the number of
    integers in [0, ..., ab] that are not admissible.
    """
    try:
        a = int(input("Enter the positive integer a: "))
        b = int(input("Enter the positive integer b: "))

        if a <= 0 or b <= 0:
            print("Error: a and b must be positive integers.")
            return

        ab = a * b
        # The range of k is [0, 1, ..., ab], which contains ab + 1 integers.
        total_integers_in_range = ab + 1

        print(f"\nFor dimensions a = {a} and b = {b}:")

        # Case 1: If a=1 or b=1, the matrices are vectors. Any non-zero matrix has rank 1.
        # Since the basis matrices must be non-zero, all 'ab' of them have rank 1.
        # Thus, k=ab is the only admissible value.
        if a == 1 or b == 1:
            # Number of admissible values is 1 (only k=ab).
            num_admissible = 1
            num_not_admissible = total_integers_in_range - num_admissible
            
            print(f"The number of non-admissible integers is given by the equation:")
            # We print the final calculation: (total integers) - (admissible count) = result
            print(f"{total_integers_in_range} - {num_admissible} = {num_not_admissible}")

        # Case 2: If a>1 and b>1, it's a known result that k can be any integer
        # in [0, ..., ab] except for k = ab-1.
        else:
            # Number of admissible values is (ab+1) - 1 = ab.
            num_admissible = ab
            num_not_admissible = total_integers_in_range - num_admissible

            print(f"The number of non-admissible integers is given by the equation:")
            # We print the final calculation: (total integers) - (admissible count) = result
            print(f"{total_integers_in_range} - {num_admissible} = {num_not_admissible}")

    except ValueError:
        print("Invalid input. Please enter valid integers for a and b.")

# Execute the function
solve_admissible_integers()