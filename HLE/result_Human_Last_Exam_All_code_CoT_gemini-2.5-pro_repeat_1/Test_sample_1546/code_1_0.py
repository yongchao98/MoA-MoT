def solve():
    """
    This function finds and prints the irreducible characters of degree 4 for S5.
    
    The symmetric group S5 has 7 conjugacy classes, and thus 7 irreducible characters.
    The degrees of these characters correspond to the partitions of 5.
    By applying the hook-length formula, we find that the partitions corresponding
    to characters of degree 4 are (4,1) and (2,1,1,1).

    The conjugacy classes of S5 are represented by permutations with cycle structures:
    1. (1,1,1,1,1) - e
    2. (2,1,1,1) - e.g., (1 2)
    3. (3,1,1) - e.g., (1 2 3)
    4. (4,1) - e.g., (1 2 3 4)
    5. (5) - e.g., (1 2 3 4 5)
    6. (2,2,1) - e.g., (1 2)(3 4)
    7. (3,2) - e.g., (1 2 3)(4 5)

    The character chi_(4,1) is the standard character of S5. Its values on the
    conjugacy classes are [4, 2, 1, 0, -1, 0, -1].

    The character chi_(2,1,1,1) is the product of chi_(4,1) and the sign character.
    Its values on the conjugacy classes are [4, -2, 1, 0, -1, 0, 1].

    The problem asks for the list of values for each character, sorted in
    ascending order.
    """

    # Values for the character corresponding to the partition (4,1)
    char1_values = [4, 2, 1, 0, -1, 0, -1]
    
    # Values for the character corresponding to the partition (2,1,1,1)
    char2_values = [4, -2, 1, 0, -1, 0, 1]

    # Sort the values in ascending order as requested
    sorted_char1 = sorted(char1_values)
    sorted_char2 = sorted(char2_values)
    
    # Print the two lists, separated by a comma
    print(f"{sorted_char1}, {sorted_char2}")

solve()