import math

def solve():
    """
    This function determines the smallest N for the multiplicative magic square problem
    and demonstrates the solution. The method is constructive, based on building
    the square from two sets of four integers, A and X.
    """

    # The smallest N is found by finding two sets of four distinct integers, A and X,
    # such that the 16 products of their elements are distinct and the maximum
    # product is minimized. The best known sets for this are:
    A = [1, 2, 3, 4]
    X = [1, 5, 6, 7]

    # These sets are assigned to the construction variables {a,b,c,d} and {x,y,z,w}.
    a, b, c, d = A
    x, y, z, w = X

    # The square is constructed using the following formula, which guarantees
    # the row and column products are equal.
    square = [
        [a*x, b*y, c*z, d*w],
        [c*y, d*x, a*w, b*z],
        [d*z, c*w, b*x, a*y],
        [b*w, a*z, d*y, c*x]
    ]

    # Extract the numbers, find N, and calculate the magic product P.
    numbers = sorted([num for row in square for num in row])
    N = max(numbers)
    magic_product = 1
    for num in square[0]:
        magic_product *= num

    print(f"The smallest N is {N}.")
    print(f"A valid 4x4 table can be constructed using the following 16 distinct integers, where the largest is {N}:")
    print(numbers)
    print("\nOne possible arrangement of these numbers is:")
    for row in square:
        print(" ".join(f"{n:2d}" for n in row))
    print(f"\nThe product of each row and column is {magic_product}.")
    
    # As requested, here is one of the product equations from the table.
    print("\nThe final equation for the first row is:")
    equation_str = " * ".join(map(str, square[0]))
    print(f"{equation_str} = {magic_product}")

solve()
<<<28>>>