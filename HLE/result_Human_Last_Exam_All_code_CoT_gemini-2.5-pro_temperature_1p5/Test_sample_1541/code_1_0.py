def solve():
    """
    Reads a sum string, sorts the numbers, and prints the new sum.
    Example: 3+2+1 -> 1+2+3
    """
    # Read the input string, e.g., "3+2+1"
    s = input()
    
    # Split the string by '+' to get a list of number strings: ['3', '2', '1']
    numbers = s.split('+')
    
    # Sort the list of numbers. Lexicographical sort works since they are single digits.
    # ['3', '2', '1'] becomes ['1', '2', '3']
    numbers.sort()
    
    # Join the sorted list back into a string with '+': "1+2+3"
    new_sum = "+".join(numbers)
    
    # Print the final result
    print(new_sum)

solve()