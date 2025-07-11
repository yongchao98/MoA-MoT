import math

def calculate_circulons():
    """
    Calculates the number of circulon defects for a gauge theory with G=SO(3)
    in (d+1) dimensions for d=1 to 6.

    The number of circulons is given by |pi_d(SO(3))| * |pi_1(SO(3))|.
    """

    # Homotopy groups of SO(3) and their sizes.
    # For d>=2, pi_d(SO(3)) = pi_d(S^3).
    # pi_d_groups maps dimension d to the string representation of the group pi_d(SO(3)).
    # pi_d_sizes maps dimension d to the size of the group |pi_d(SO(3))|.
    pi_d_groups = {
        1: "Z_2",
        2: "0",
        3: "Z",
        4: "Z_2",
        5: "Z_2",
        6: "Z_12",
    }
    pi_d_sizes = {
        1: 2,
        2: 1,  # The trivial group {0} has size 1
        3: math.inf, # The integers Z are an infinite group
        4: 2,
        5: 2,
        6: 12,
    }

    # The fundamental group pi_1(SO(3)) is Z_2.
    pi_1_group_str = "Z_2"
    pi_1_size = 2

    print("Calculating the number of circulons for G=SO(3) in d+1 dimensions.")
    print("Formula: Number = |pi_d(SO(3))| * |pi_1(SO(3))|\n")

    for d in range(1, 7):
        group_d_str = pi_d_groups[d]
        size_d = pi_d_sizes[d]

        if size_d == math.inf:
            result_str = "infinity"
            size_d_str = "infinity"
        else:
            result = int(size_d * pi_1_size)
            result_str = str(result)
            size_d_str = str(int(size_d))
        
        # Print the full equation for each dimension d
        print(f"For d={d}:")
        print(f"Number = |pi_{d}(SO(3))| * |pi_1(SO(3))|")
        print(f"       = |{group_d_str}| * |{pi_1_group_str}|")
        print(f"       = {size_d_str} * {pi_1_size} = {result_str}\n")

calculate_circulons()