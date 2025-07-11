def solve():
    """
    Finds all integers n between 7 and 55 (inclusive) such that n is not
    divisible by 7.
    """
    possible_n_values = []
    for n in range(7, 56):
        if n % 7 != 0:
            possible_n_values.append(n)
    
    # The prompt asks to output each number in the final equation.
    # We will interpret this as printing the list of numbers.
    print(*possible_n_values)

solve()