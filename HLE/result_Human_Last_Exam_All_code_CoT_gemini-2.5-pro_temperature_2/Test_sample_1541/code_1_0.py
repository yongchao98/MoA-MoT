import sys

def solve():
    """
    This function reads a sum as a string, rearranges the numbers in non-decreasing
    order, and prints the new sum.
    Example: 1+3+2+1 -> 1+1+2+3
    """
    # Read the sum from standard input.
    s = sys.stdin.readline().strip()

    # The summands are separated by '+'. We can split the string by the '+'
    # character to get a list of the numbers (as strings).
    numbers = s.split('+')

    # Sort the list of numbers. Since the numbers are single digits '1', '2', '3', '4',
    # a simple lexicographical sort works correctly.
    numbers.sort()

    # Join the sorted list of numbers back into a single string, with '+'
    # as the separator.
    new_sum = '+'.join(numbers)

    # Print the final rearranged sum.
    print(new_sum)

solve()