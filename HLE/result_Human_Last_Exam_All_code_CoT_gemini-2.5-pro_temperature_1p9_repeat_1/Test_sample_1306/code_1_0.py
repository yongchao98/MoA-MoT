import sys

def solve():
    """
    This function calculates the percentage of irreducible representations among the
    indecomposable tilting modules for the quantum group u_q(sl_2) where q is
    a primitive 3rd root of unity.
    """

    # For q a primitive l-th root of unity, l=3.
    l = 3

    # The number of irreducible representations (which are also tilting modules) is l - 1.
    num_irreducible = l - 1

    # The total number of indecomposable tilting modules is 2 * (l - 1),
    # consisting of the l-1 simple modules and their l-1 projective covers.
    total_indecomposable_tilting = 2 * (l - 1)

    # Calculate the percentage.
    percentage = (num_irreducible / total_indecomposable_tilting) * 100

    # Print the equation as requested.
    # The numbers in the final equation are num_irreducible and total_indecomposable_tilting.
    print(f"Under the interpretation of considering the category of tilting modules:")
    print(f"Number of irreducible modules: {num_irreducible}")
    print(f"Total number of indecomposable tilting modules: {total_indecomposable_tilting}")
    print(f"The percentage is ({num_irreducible} / {total_indecomposable_tilting}) * 100 = {int(percentage)}%")


solve()