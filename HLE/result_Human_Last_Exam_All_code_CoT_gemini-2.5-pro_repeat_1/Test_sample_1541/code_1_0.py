import sys

def solve():
    """
    Reads a sum string, sorts the numbers, and prints the new sum.
    Example: 3+2+1 -> 1+2+3
    """
    # Read the input string from the user.
    # The .strip() removes any potential leading/trailing whitespace, including the newline.
    s = sys.stdin.readline().strip()

    # A memory-efficient way is to count the occurrences of each number.
    count1 = s.count('1')
    count2 = s.count('2')
    count3 = s.count('3')
    count4 = s.count('4') # As per the prompt's problem description

    # Create a list that will hold the numbers in sorted order.
    sorted_numbers = []

    # Extend the list with the correct number of '1's, '2's, '3's, and '4's.
    sorted_numbers.extend(['1'] * count1)
    sorted_numbers.extend(['2'] * count2)
    sorted_numbers.extend(['3'] * count3)
    sorted_numbers.extend(['4'] * count4)

    # Join the sorted numbers with '+' to form the final sum string.
    # The output will correctly show each number in the final equation.
    # For example, ['1', '1', '2'] becomes "1+1+2".
    result = "+".join(sorted_numbers)

    # Print the final result.
    print(result)

solve()