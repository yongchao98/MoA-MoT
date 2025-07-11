def find_winning_t_for_p2(limit):
    """
    This function finds the values of T for which the second player has a winning strategy
    in the described token game. These values are the Fibonacci numbers.

    Args:
    limit (int): The maximum value for T to check.

    Returns:
    list: A list of Fibonacci numbers up to the given limit.
    """
    # The set of losing positions for the first player are the Fibonacci numbers {1, 2, 3, 5, ...}.
    # We can generate this sequence starting with a=1 and b=2.
    winning_values = []
    a, b = 1, 2
    
    # Handle the first value T=1
    if 1 <= limit:
        winning_values.append(a)
    else:
        return []

    # Generate subsequent Fibonacci numbers
    while b <= limit:
        winning_values.append(b)
        a, b = b, a + b
        
    return winning_values

if __name__ == '__main__':
    # Set an upper bound for T to demonstrate the pattern.
    T_limit = 1000
    
    winning_T = find_winning_t_for_p2(T_limit)
    
    print("The second player has a winning strategy if and only if the initial number of tokens, T, is a Fibonacci number.")
    print(f"The winning values for T up to {T_limit} are:")
    
    # The prompt requires printing each number in the final equation. As there is no equation,
    # we will just print the resulting numbers.
    # To format the output nicely, we can convert the list to a string.
    winning_T_str = ', '.join(map(str, winning_T))
    print(winning_T_str)