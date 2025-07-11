import math

def solve_wuxing_factorial():
    """
    Calculates the memory cost (z) for an optimized Wuxing C program
    to compute 100! and finds the first 3 digits (y) of the result.
    """

    # --- Step 1: Analyze the result of 100! ---

    # Calculate 100!
    factorial_100 = math.factorial(100)
    
    # Get the number of decimal digits in 100!
    num_digits = len(str(factorial_100))
    
    # Get the first 3 digits (y)
    y = int(str(factorial_100)[:3])

    print("Analyzing the problem:")
    print(f"100! has {num_digits} decimal digits.")
    print(f"The first 3 digits ('y') of 100! are: {y}\n")
    

    # --- Step 2: Calculate the optimized memory size (z) ---
    
    print("Calculating minimum memory usage ('z') for variables in an optimized Wuxing C program:")
    
    # Define memory sizes for XVM data types in Decimal Digits (D)
    D_digit = 1
    D_cent = 2
    D_char = 3
    
    # To store the 158-digit result, the most memory-efficient array is one
    # where each element stores a single digit.
    # Variable: digit result[158];
    mem_result_array = num_digits * D_digit
    print(f"- Memory for result array: A {num_digits}-element array of 'digit' (1D) is needed. Size = {num_digits} * {D_digit}D = {mem_result_array}D")

    # The main loop counter 'i' goes from 2 to 100.
    # The smallest type that can hold 100 is 'cent' (0-99) or 'char' (0-999). 
    # Since 100 > 99, we must use 'char'. A loop to 100 implies `i <= 100` so 100 must be storable.
    # Oh wait, a standard 'for' loop condition `for(i=2; i<=100; i++)` uses the type 'cent' correctly.
    # When `i` is 99, it is incremented to 100, the comparison `100 <= 100` is true. `cent` is sufficient.
    # Re-evaluating: max value of i is 100. 'cent' range is 0-99. It cannot hold 100.
    # The next size up is 'char' (3D), which holds 0-999.
    # Variable: char i;
    mem_i = D_char
    print(f"- Memory for loop counter 'i' (up to 100): Needs to hold 100. 'cent' (max 99) is too small. 'char' (3D) is required. Size = {mem_i}D")

    # Inner loop counter 'j' iterates through the digits of the result array.
    # It goes up to `num_digits - 1`, so max value is 157.
    # 'char' (0-999) is required.
    # Variable: char j;
    mem_j = D_char
    print(f"- Memory for inner loop counter 'j' (up to {num_digits-1}): 'char' (3D) is required. Size = {mem_j}D")
    
    # A variable is needed to track the current number of digits in the result.
    # Max value is 158.
    # Variable: char max_size;
    mem_max_size = D_char
    print(f"- Memory for 'max_size' tracker (up to {num_digits}): 'char' (3D) is required. Size = {mem_max_size}D")

    # In the multiplication `product = result[j] * i + carry`, we need to find the
    # max value of 'carry'. The multiplication is by a number 'i' up to 100.
    # For a digit d, the operation is d * i + carry.
    # Max product_step = 9 * 100 + previous_carry.
    # Max previous_carry = floor((9 * 100 + 99)/10) = 99.
    # So max carry is 99. 'cent' (0-99) is sufficient.
    # Variable: cent carry;
    mem_carry = D_cent
    print(f"- Memory for 'carry' (max value is 99): 'cent' (2D) is sufficient. Size = {mem_carry}D")

    # The temporary product `temp = result[j] * i + carry`.
    # Max value = 9 * 100 + 99 = 999.
    # 'char' (0-999) is required and perfectly fits.
    # Variable: char temp;
    mem_temp = D_char
    print(f"- Memory for temporary 'temp' product (max value is 999): 'char' (3D) is sufficient. Size = {mem_temp}D")

    # Sum all memory costs for 'z'
    z = mem_result_array + mem_i + mem_j + mem_max_size + mem_carry + mem_temp

    print("\n--- Final Equation ---")
    print(f"z = (mem_result) + (mem_i) + (mem_j) + (mem_max_size) + (mem_carry) + (mem_temp)")
    print(f"z = {mem_result_array}D + {mem_i}D + {mem_j}D + {mem_max_size}D + {mem_carry}D + {mem_temp}D = {z}D")

    print("\n--- Result ---")
    print(f"{z}:{y}")


solve_wuxing_factorial()