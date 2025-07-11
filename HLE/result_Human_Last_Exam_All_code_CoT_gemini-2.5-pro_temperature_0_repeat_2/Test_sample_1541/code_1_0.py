def solve():
    """
    Reads a sum of numbers (1, 2, 3, 4), sorts the numbers,
    and prints the new sum.
    """
    # Read the input string, e.g., "3+2+1"
    s = input()

    # Split the string by the '+' character to get a list of numbers as strings
    # e.g., ["3", "2", "1"]
    numbers = s.split('+')

    # Sort the list of numbers. Since they are single digits,
    # lexical sorting is equivalent to numerical sorting.
    # e.g., ["1", "2", "3"]
    numbers.sort()

    # Join the sorted list of numbers back into a string with '+' as the separator
    # e.g., "1+2+3"
    new_sum = "+".join(numbers)

    # Print the final rearranged sum
    print(new_sum)

solve()