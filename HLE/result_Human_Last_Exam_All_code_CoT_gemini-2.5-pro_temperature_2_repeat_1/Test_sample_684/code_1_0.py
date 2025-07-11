def solve():
    """
    Finds and prints the values of n in the range [7, 55] for which it is
    possible to leave just one gift on the 5D hypercube.

    The problem can be solved by analyzing the underlying linear algebra over F_2.
    A configuration is solvable if and only if n satisfies the condition n % 7 == 1.
    """

    possible_n = []
    for n in range(7, 56):
        # The condition for solvability is that n mod 7 equals 1.
        if n % 7 == 1:
            possible_n.append(n)

    # Print the numbers in increasing order
    print("The possible values for n are:")
    for n in possible_n:
        print(n)

solve()