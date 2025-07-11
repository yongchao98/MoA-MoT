def solve():
    """
    Finds and prints all values of n in the range [7, 55] for which it is
    possible to reach a state with a single gift.
    The condition for possibility is that n is not a multiple of 7.
    """
    possible_n = []
    for n in range(7, 56):
        if n % 7 != 0:
            possible_n.append(n)
    
    # The problem asks for the values in increasing order.
    # The loop naturally finds them in increasing order.
    # To match the requested format, we convert them to strings for printing.
    result_str = ", ".join(map(str, possible_n))
    print(result_str)
    print(f"<<<{result_str}>>>")

solve()