def solve_admissible_matrices(a, b):
    """
    Calculates the number of non-admissible integers k for given matrix dimensions a and b.

    An integer k is "admissible" if there exist ab non-zero complex a x b matrices A_i
    that are mutually orthogonal (tr(A_i^dagger A_j)=0) and exactly k of which have rank 1.
    This function finds the number of integers in {0, 1, ..., ab} that are not admissible.
    """
    if not (isinstance(a, int) and isinstance(b, int) and a > 0 and b > 0):
        print("Please provide positive integers for a and b.")
        return

    print(f"Solving for a={a}, b={b}:")

    # Case 1: The matrices are essentially vectors.
    if min(a, b) == 1:
        # In this case, any non-zero a x b matrix must have rank 1.
        # The problem requires all ab matrices in the basis to be non-zero,
        # so all of them must have rank 1.
        # This means the only admissible value for k is a*b.
        # The non-admissible integers are {0, 1, ..., a*b - 1}.
        result = a * b
        
        print("This is the case where min(a, b) = 1.")
        print("In this scenario, all non-zero matrices have rank 1.")
        print(f"Therefore, the only admissible k is k = a*b = {result}.")
        print(f"The number of non-admissible integers is the count of numbers from 0 to {result-1}.")
        print(f"Final calculation: a * b = {a} * {b} = {result}")

    # Case 2: The general case for matrices with dimensions >= 2.
    else:
        prod = a * b
        print(f"This is the case where min(a, b) >= 2. The product a*b is {prod}.")
        print("The number of rank-1 matrices in any basis must have the same parity as a*b.")
        print("The non-admissible integers are therefore those with a different parity than a*b.")

        # If a*b is even, the non-admissible integers are the odd ones.
        if prod % 2 == 0:
            result = prod // 2
            print(f"Since a*b is even, the non-admissible integers are the odd numbers in the range 0...{prod}.")
            print(f"Final calculation: (a*b) / 2 = {prod} / 2 = {result}")
            
        # If a*b is odd, the non-admissible integers are the even ones.
        else:
            result = (prod + 1) // 2
            print(f"Since a*b is odd, the non-admissible integers are the even numbers in the range 0...{prod}.")
            print(f"Final calculation: (a*b + 1) / 2 = ({prod} + 1) / 2 = {result}")

    print("-" * 20)


# You can change the values of a and b here to test different cases.
a = 4
b = 5
solve_admissible_matrices(a, b)

a = 3
b = 5
solve_admissible_matrices(a, b)

a = 1
b = 7
solve_admissible_matrices(a, b)