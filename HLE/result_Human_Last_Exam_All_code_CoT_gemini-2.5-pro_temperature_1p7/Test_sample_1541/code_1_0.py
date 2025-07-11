def solve():
    """
    Reads a sum of numbers (1, 2, 3, 4), sorts the numbers,
    and prints the new sum in non-decreasing order.
    """
    # 1. Read the input sum from a single line
    s = input()

    # 2. Extract the numbers by splitting the string at the '+' signs
    summands = s.split('+')

    # 3. Sort the list of numbers. Since they are single digits,
    # alphabetical sorting is equivalent to numerical sorting.
    summands.sort()

    # 4. Join the sorted numbers back into a string with '+'
    new_sum = "+".join(summands)

    # 5. Print the final result. The output format is the full equation.
    print(new_sum)

solve()