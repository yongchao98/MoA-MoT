import sys

def solve():
    """
    Reads a sum, sorts the numbers, and prints the new sum.
    """
    # Read the input line, e.g., "3+2+1"
    s = sys.stdin.readline().strip()
    
    # Split the string by '+' to get a list of numbers, e.g., ['3', '2', '1']
    numbers = s.split('+')
    
    # Sort the list of numbers. For single digits, default sort is fine.
    # e.g., ['1', '2', '3']
    numbers.sort()
    
    # Join the sorted numbers back with '+' and print
    # e.g., "1+2+3"
    print('+'.join(numbers))

solve()