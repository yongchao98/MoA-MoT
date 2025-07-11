def solve():
    """
    Reads a sum as a string, sorts the numbers, and prints the new sum.
    """
    # Read the input string, e.g., "3+2+1+4"
    s = input()

    # The numbers are separated by '+'. We can split the string to get them.
    # Since '1' < '2' < '3' < '4' lexicographically, we don't need to convert to integers.
    summands = s.split('+')

    # Sort the list of numbers.
    summands.sort()

    # Join the sorted list back into a string with '+' as the separator.
    # The final equation will have each number printed.
    new_sum = '+'.join(summands)

    # Print the final result.
    print(new_sum)

solve()