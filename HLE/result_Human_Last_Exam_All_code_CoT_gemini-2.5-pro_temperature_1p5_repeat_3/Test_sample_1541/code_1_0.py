def solve():
    """
    Reads a sum of numbers (1, 2, 3, 4), rearranges them in non-decreasing order,
    and prints the new sum.
    """
    # Read the input string
    s = input()

    # Extract the numbers by splitting the string at the '+' character
    numbers_str = s.split('+')

    # Sort the list of numbers. Python's sort on strings of single digits works
    # correctly, so we don't strictly need to convert to integers, but doing
    # so is more robust if numbers could have more than one digit.
    numbers_str.sort()

    # Join the sorted numbers back into a single string with '+'
    result = '+'.join(numbers_str)

    # Print the final rearranged sum
    print(result)

solve()