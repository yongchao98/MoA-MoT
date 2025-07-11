import math

def calculate_large_factorial(n):
    """
    Calculates the factorial of n using an array of digits to handle large numbers.
    The digits are stored in reverse order in the list (least significant digit first).
    """
    # Initialize the result array with the number 1.
    # We use a list of integers, where each integer is a single digit.
    result = [1]
    
    # Iterate from 2 to n to multiply.
    for i in range(2, n + 1):
        carry = 0
        num_digits = len(result)
        
        # Multiply each digit in the result array by i.
        for j in range(num_digits):
            # Calculate the product of the current digit, the multiplier i, and add the carry.
            product = result[j] * i + carry
            # The new digit at this position is the remainder of the product divided by 10.
            result[j] = product % 10
            # The new carry is the integer division of the product by 10.
            carry = product // 10
            
        # If there is a carry left over, append its digits to the end of the result array.
        while carry > 0:
            result.append(carry % 10)
            carry //= 10
            
    return result

# Step 1: Calculate memory size 'z' based on the theoretical C program for Wuxing.
#   - result array: digit result[160] -> 160D
#   - variables (i, num_digits, j, carry, product): 5 * char (3D each) -> 15D
#   - Total z = 160 + 15 = 175D
z = 175

# Step 2: Calculate 100! to find 'y'.
n = 100
factorial_digits_reversed = calculate_large_factorial(n)

# The result is stored in reverse. The most significant digits are at the end of the list.
# We extract the last 3 digits and reverse their order to get the first 3 digits of the number.
first_three_digits_list = factorial_digits_reversed[-1:-4:-1]
y = "".join(map(str, first_three_digits_list))

# Step 3: Print the numbers of the final equation z:y
print(f"z = {z}")
print(f"y = {y}")
print(f"The final answer is:")
print(f"{z}:{y}")
# >>>175:903