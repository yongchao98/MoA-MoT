import sys

def solve():
    """
    This function reads a sum from standard input, sorts the numbers,
    and prints the new sum in non-decreasing order.
    Example: 1+3+2+1 -> 1+1+2+3
    """
    # Read the entire line of input, e.g., "1+3+2+1"
    s = sys.stdin.readline().strip()

    # The input is guaranteed to be non-empty.
    # Split the string by the '+' character to get a list of number strings.
    # e.g., "1+3+2+1" becomes ['1', '3', '2', '1']
    numbers = s.split('+')

    # Sort the list of number strings. Python's sort for strings works
    # correctly here since they are single digits.
    # e.g., ['1', '3', '2', '1'] becomes ['1', '1', '2', '3']
    numbers.sort()

    # Join the sorted list of strings back together with '+' as the separator.
    # e.g., ['1', '1', '2', '3'] becomes "1+1+2+3"
    new_sum = "+".join(numbers)

    # Print the final result.
    print(new_sum)

solve()