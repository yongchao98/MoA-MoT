def solve_xenia_sum():
    """
    Reads a sum of numbers (1, 2, 3, 4), sorts the numbers,
    and prints the new sum.
    """
    # Read the input string, e.g., "3+2+1"
    s = input()

    # Split the string by the '+' delimiter to get a list of number strings
    # e.g., "3+2+1" becomes ['3', '2', '1']
    numbers = s.split('+')

    # Sort the list of numbers. Since they are single digits,
    # lexicographical sorting is equivalent to numerical sorting.
    # e.g., ['3', '2', '1'] becomes ['1', '2', '3']
    numbers.sort()

    # Join the sorted list back into a single string with '+'
    # e.g., ['1', '2', '3'] becomes "1+2+3"
    sorted_sum = "+".join(numbers)

    # Print the final rearranged sum
    print(sorted_sum)

solve_xenia_sum()