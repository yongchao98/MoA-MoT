def solve():
    """
    Reads a sum of numbers (1, 2, 3, 4), rearranges them in non-decreasing order,
    and prints the new sum.
    """
    s = input()

    # The numbers are only 1, 2, 3, or 4. We can use counting sort.
    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0

    # Iterate through the string, stepping by 2 to only get the digits.
    for i in range(0, len(s), 2):
        char = s[i]
        if char == '1':
            count1 += 1
        elif char == '2':
            count2 += 1
        elif char == '3':
            count3 += 1
        elif char == '4':
            count4 += 1

    # Build the result list of strings in the correct order.
    result_parts = []
    result_parts.extend(['1'] * count1)
    result_parts.extend(['2'] * count2)
    result_parts.extend(['3'] * count3)
    result_parts.extend(['4'] * count4)

    # Join the parts with '+' and print the final equation.
    # The problem requires outputting each number in the final equation.
    # The join function handles this by placing '+' between each number.
    print("+".join(result_parts))

solve()