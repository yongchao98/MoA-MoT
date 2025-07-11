import math

def get_character_values():
    """
    Calculates the two irreducible characters of degree 4 for the symmetric group S5.
    """
    # The 7 conjugacy classes of S5, represented by their cycle structures (partitions of 5)
    class_partitions = [
        [1, 1, 1, 1, 1],  # Identity (e.g., ())
        [2, 1, 1, 1],      # Transposition (e.g., (1 2))
        [3, 1, 1],        # 3-cycle (e.g., (1 2 3))
        [4, 1],           # 4-cycle (e.g., (1 2 3 4))
        [5],              # 5-cycle (e.g., (1 2 3 4 5))
        [2, 2, 1],        # Double transposition (e.g., (1 2)(3 4))
        [3, 2]            # 3-cycle and 2-cycle (e.g., (1 2 3)(4 5))
    ]

    # --- Character 1: Corresponding to partition [4,1] ---
    # This character is given by chi_standard - chi_trivial.
    # chi_standard(g) = number of fixed points of permutation g.
    # chi_trivial(g) = 1 for all g.
    char_values_1 = []
    for p in class_partitions:
        # The number of fixed points is the number of 1s in the cycle partition.
        fixed_points = p.count(1)
        value = fixed_points - 1
        char_values_1.append(value)

    # --- Character 2: Corresponding to partition [2,1,1,1] ---
    # This character is (chi_[4,1]) x (sign_character).
    # The sign of a permutation with n elements and k cycles is (-1)^(n - k).
    char_values_2 = []
    n = 5
    for i, p in enumerate(class_partitions):
        num_cycles = len(p)
        # Permutation is even if (n - num_cycles) is even, odd otherwise.
        sign = (-1)**(n - num_cycles)
        value = char_values_1[i] * sign
        char_values_2.append(value)

    # Sort the character values in ascending order as requested
    sorted_char_1 = sorted(char_values_1)
    sorted_char_2 = sorted(char_values_2)
    
    # Print the final result in the specified format
    print(f"{sorted_char_1}, {sorted_char_2}")

get_character_values()