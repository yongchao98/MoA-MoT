import math

def solve_circulons():
    """
    Calculates the number of circulon types for dimensions d=1 to d=6.

    A circulon is a circular defect in a d+1 dimensional SO(3) gauge theory.
    Its classification is given by the homotopy classes of maps [S^1 x S^(d-1), SO(3)].
    The number of classes is |pi_1(SO(3))| * |pi_{d-1}(SO(3))| * |pi_d(SO(3))|.
    """

    # We store the size of the homotopy groups of SO(3).
    # -1 represents an infinite group (like Z).
    # For k > 1, pi_k(SO(3)) = pi_k(S^3).
    pi_sizes = {
        0: 1,    # |pi_0(SO(3))|, from connectedness
        1: 2,    # |pi_1(SO(3))| = |Z_2|
        2: 1,    # |pi_2(SO(3))| = |pi_2(S^3)| = 0
        3: -1,   # |pi_3(SO(3))| = |pi_3(S^3)| = |Z|
        4: 2,    # |pi_4(SO(3))| = |pi_4(S^3)| = |Z_2|
        5: 2,    # |pi_5(SO(3))| = |pi_5(S^3)| = |Z_2|
        6: 12,   # |pi_6(SO(3))| = |pi_6(S^3)| = |Z_12|
    }

    # Helper function to format the size for printing
    def size_to_str(size):
        return "infinity" if size == -1 else str(size)

    print("The number of circulon types, N(d), is calculated for d=1 to 6.")
    print("Formula: N(d) = |π₁(SO(3))| * |πₑ₋₁(SO(3))| * |πₑ(SO(3))|\n")

    for d in range(1, 7):
        p1_size = pi_sizes[1]
        pd_minus_1_size = pi_sizes[d - 1]
        pd_size = pi_sizes[d]
        
        # Calculate the result
        if -1 in [p1_size, pd_minus_1_size, pd_size]:
            result_str = "infinity"
        else:
            result = p1_size * pd_minus_1_size * pd_size
            result_str = str(result)
            
        # Print the equation with numbers substituted
        term1 = size_to_str(p1_size)
        term2 = size_to_str(pd_minus_1_size)
        term3 = size_to_str(pd_size)
        
        print(f"For d={d}:")
        print(f"N({d}) = |π₁| * |π_{d-1}| * |π_{d}| = {term1} * {term2} * {term3} = {result_str}")

# Run the calculation
solve_circulons()