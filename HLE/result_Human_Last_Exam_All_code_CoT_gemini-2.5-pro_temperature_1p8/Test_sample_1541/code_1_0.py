import sys

def solve():
    """
    Reads a string representing a sum of numbers (1, 2, 3, 4),
    sorts the numbers, and prints the new sum.
    """
    # Read the sum from standard input
    s = sys.stdin.readline().strip()
    
    # Split the string by the '+' delimiter to get a list of the numbers
    summands = s.split('+')
    
    # Sort the list of numbers. Since they are single digits,
    # default string sort works correctly.
    summands.sort()
    
    # Join the sorted numbers back together with '+'
    result = "+".join(summands)
    
    # Print the final rearranged sum
    print(result)

solve()