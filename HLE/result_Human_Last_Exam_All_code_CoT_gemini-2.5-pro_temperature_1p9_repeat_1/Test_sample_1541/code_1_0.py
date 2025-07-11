def solve():
    """
    Reads a string representing a sum of numbers 1, 2, 3, and 4,
    rearranges the numbers in non-decreasing order, and prints the new sum.
    """
    # Read the input string, for example "3+2+1"
    s = input()

    # Split the string by the '+' character to get a list of the numbers.
    # For "3+2+1", this results in ['3', '2', '1']
    numbers = s.split('+')

    # Sort the list of numbers.
    # For ['3', '2', '1'], this becomes ['1', '2', '3']
    numbers.sort()

    # Join the sorted list back into a string with '+' as the separator.
    # For ['1', '2', '3'], this becomes "1+2+3"
    output_sum = "+".join(numbers)

    # Print the final result. The output needs to be the full equation.
    print(output_sum)

solve()