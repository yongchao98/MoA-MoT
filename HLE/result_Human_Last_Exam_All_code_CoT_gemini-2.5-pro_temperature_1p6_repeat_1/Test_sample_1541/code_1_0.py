def solve():
    """
    Reads a sum as a string, sorts the numbers, and prints the new sum.
    """
    # Read the input string, e.g., "3+2+1"
    s = input()

    # Extract all the numbers from the string.
    # The numbers are the characters that are not '+'
    numbers = []
    for char in s:
        if char != '+':
            numbers.append(char)

    # Sort the list of numbers. Since they are single digits,
    # sorting them as strings works correctly.
    numbers.sort()

    # Join the sorted numbers with '+' and print the result.
    # The problem asks to print the full equation.
    output = ""
    for i in range(len(numbers)):
        output += numbers[i]
        if i < len(numbers) - 1:
            output += "+"
    print(output)

solve()