def solve():
    """
    Finds and prints all values of n between 7 and 55 (inclusive) for which
    it is possible to leave a single gift on a 5-D hypercube.
    
    The condition for possibility is that n is not divisible by 7.
    """
    
    possible_n_values = []
    for n in range(7, 56):
        if n % 7 != 0:
            possible_n_values.append(n)
            
    print(*possible_n_values)

solve()