import sys

def solve():
    """
    Reads a sum like "1+3+2", sorts the numbers, and prints the new sum.
    """
    # Read the sum from standard input
    s = sys.stdin.readline().strip()
    
    # The sum contains numbers separated by '+'. We extract them.
    summands = s.split('+')
    
    # Sort the numbers in non-decreasing order.
    # Since they are single digits, lexicographical sort is equivalent to numerical sort.
    summands.sort()
    
    # Join the sorted numbers back together with '+' to form the new sum.
    new_sum = "+".join(summands)
    
    # Print the final rearranged sum. The problem states to output each number
    # in the final equation.
    print(new_sum)

solve()