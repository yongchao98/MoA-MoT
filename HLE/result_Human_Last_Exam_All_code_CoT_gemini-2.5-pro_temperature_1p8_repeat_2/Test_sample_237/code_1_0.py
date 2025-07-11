def calculate_factorial_and_memory():
    """
    Calculates the value of 100! using large number multiplication
    and determines the optimized memory usage based on the Wuxing architecture.
    """

    # Part 1: Memory size (z) calculation, as determined by the plan.
    # result[160] * 1D = 160D
    # i * 3D (char)   = 3D
    # j * 3D (char)   = 3D
    # result_size * 3D (char) = 3D
    # carry * 3D (char) = 3D
    # product * 3D (char) = 3D
    # Total = 160 + 3 + 3 + 3 + 3 + 3 = 175
    z = 175

    # Part 2: Factorial calculation to find y
    
    # Initialize result list to represent the number 1
    result_digits = [1]
    
    # Multiply by numbers from 2 to 100
    for i in range(2, 101):
        carry = 0
        # Iterate through each digit of the current result
        for j in range(len(result_digits)):
            product = result_digits[j] * i + carry
            result_digits[j] = product % 10
            carry = product // 10
            
        # If there's a remaining carry, append its digits to the list
        while carry > 0:
            result_digits.append(carry % 10)
            carry = carry // 10

    # The result_digits list stores the number in reverse order.
    # Get the first 3 digits from the end of the list.
    d1 = result_digits[-1]
    d2 = result_digits[-2]
    d3 = result_digits[-3]
    
    y = f"{d1}{d2}{d3}"

    # Print the final answer in the format z:y
    print(f"{z}:{y}")

# Execute the function
calculate_factorial_and_memory()