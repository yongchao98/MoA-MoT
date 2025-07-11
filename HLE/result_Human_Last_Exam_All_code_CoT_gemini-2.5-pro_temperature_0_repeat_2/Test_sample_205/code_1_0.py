def calculate_large_factorial():
    """
    This function calculates 100! using an array-based multiplication
    algorithm, similar to how it would be implemented on a memory-constrained
    system like the Wuxing computer.
    """
    n = 100
    # Initialize the result array with the starting value 1
    result = [1]

    # Loop from 2 to n, multiplying the result by each number
    for i in range(2, n + 1):
        carry = 0
        # Iterate through each digit of the current result
        for j in range(len(result)):
            # Perform multiplication and add the carry from the previous step
            product = result[j] * i + carry
            # The new digit is the remainder of the product
            result[j] = product % 10
            # The new carry is the integer division of the product
            carry = product // 10

        # If a carry remains, append its digits to the result array
        while carry > 0:
            result.append(carry % 10)
            carry = carry // 10

    # The result is stored in reverse order, so we print it reversed.
    # This fulfills the requirement to output each number in the final equation.
    print("100! = ", end="")
    for digit in reversed(result):
        print(digit, end="")
    print()

calculate_large_factorial()