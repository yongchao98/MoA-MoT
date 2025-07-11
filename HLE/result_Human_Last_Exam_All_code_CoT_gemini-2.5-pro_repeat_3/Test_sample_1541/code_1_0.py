def solve():
    """
    Reads a sum of numbers 1, 2, 3, 4, sorts the numbers,
    and prints the new sum.
    """
    # Read the input sum as a string
    s = input()

    # The numbers are separated by '+'. We can split the string by '+',
    # which gives us a list of the numbers as strings.
    summands = s.split('+')

    # Sort the list of numbers. Since they are all single-digit numbers,
    # lexical sorting is equivalent to numerical sorting.
    summands.sort()

    # Join the sorted numbers back together with '+' and print the result.
    # This outputs the final equation containing each sorted number.
    print('+'.join(summands))

solve()