import sys

def solve():
    """
    Reads a sum of numbers (1, 2, 3, 4), rearranges them in non-decreasing order,
    and prints the new sum.
    """
    # Read the input string, e.g., "3+2+1"
    input_sum = sys.stdin.readline().strip()
    
    # Extract the numbers from the string by filtering out the '+' characters
    numbers = [char for char in input_sum if char != '+']
    
    # Sort the list of numbers. Sorting characters '1', '2', '3' works as expected.
    numbers.sort()
    
    # Join the sorted numbers back together with '+' to form the final string
    output_sum = "+".join(numbers)
    
    # Print the result
    print(output_sum)

solve()