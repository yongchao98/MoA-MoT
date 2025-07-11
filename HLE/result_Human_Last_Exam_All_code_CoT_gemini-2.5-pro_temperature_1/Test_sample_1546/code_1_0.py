import sys

def solve():
    """
    This script calculates and prints the irreducible characters of degree 4 for the symmetric group S5.

    The two irreducible characters of degree 4 correspond to the partitions (4,1) and (2,1,1,1).

    1. Character for partition (4,1):
    The character value is (number of fixed points) - 1. The values for the 7 conjugacy
    classes of S5 are [4, 2, 1, 0, 0, -1, -1].

    2. Character for partition (2,1,1,1):
    This character is the product of the first character and the sign character. The values are
    [4, -2, 1, 0, 0, 1, -1].

    The script sorts these lists of values and prints them.
    """

    # The set of values for the character corresponding to the partition (4,1)
    char1_values = [4, 2, 1, 0, 0, -1, -1]

    # The set of values for the character corresponding to the partition (2,1,1,1)
    char2_values = [4, -2, 1, 1, 0, 0, -1]

    # Sort each list in ascending order
    char1_sorted = sorted(char1_values)
    char2_sorted = sorted(char2_values)

    # Print the final result in the specified format
    # The f-string conversion of a list automatically includes brackets and commas.
    print(f"{char1_sorted},{char2_sorted}")

solve()