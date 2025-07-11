def solve():
    """
    Reads a sum of numbers (1, 2, 3, 4), sorts the numbers,
    and prints the new sum.
    """
    # Read the input sum as a string
    s = input()

    # The summands are the characters at even indices, or more easily,
    # we can split the string by the '+' character.
    summands = s.split('+')

    # Sort the list of summands in non-decreasing order.
    # Since they are strings '1', '2', '3', '4', lexical sort is sufficient.
    summands.sort()

    # Join the sorted summands back into a single string, separated by '+'.
    # This forms the new equation that Xenia can calculate.
    new_sum = "+".join(summands)

    # Print the final rearranged sum.
    print(new_sum)

solve()