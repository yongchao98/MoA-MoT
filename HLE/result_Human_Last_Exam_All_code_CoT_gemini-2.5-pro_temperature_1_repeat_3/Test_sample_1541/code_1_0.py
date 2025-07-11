import sys

def solve():
    """
    Reads a sum of numbers (1, 2, 3, 4), sorts the numbers,
    and prints the new sum.
    """
    # Read the input sum from a single line.
    s = sys.stdin.readline().strip()

    # The numbers are separated by '+'. We can split the string
    # to get a list of the numbers (as strings).
    numbers = s.split('+')

    # Sort the list of numbers. Python's sort works lexicographically,
    # which is correct for single-digit numbers.
    numbers.sort()

    # Join the sorted numbers back into a single string,
    # separated by '+'.
    new_sum = "+".join(numbers)

    # Print the final rearranged sum.
    print(new_sum)

solve()