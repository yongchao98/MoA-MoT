import math

def solve_wuxing_factorial():
    """
    This function calculates the two parts of the problem:
    1. y: The first 3 digits of 100!
    2. z: The smallest memory size in D to calculate 100! on the Wuxing computer.
    """

    # Part 1: Calculate 'y', the first 3 digits of 100!
    # Python's math.factorial can handle large integers, making this straightforward.
    factorial_100 = math.factorial(100)
    y = str(factorial_100)[:3]

    # Part 2: Calculate 'z', the smallest memory size.
    # This is based on designing the most memory-efficient C program for the XVM.

    # First, determine the number of decimal digits in 100!. This determines
    # the size of our bignum array.
    num_digits_in_result = len(str(factorial_100)) # This is 158

    # Define the memory sizes of XVM data types in Decimal Digits (D).
    SIZES = {
        "digit": 1,
        "cent": 2,
        "char": 3,
        "int": 6,
        "long": 12
    }

    # The most memory-efficient design uses a base-10 approach with a 'digit' array.
    # We now determine the smallest possible data type for each variable required.
    
    # Variable 1: The array to hold the result digits.
    # Needs to hold 158 individual digits.
    # C declaration: digit result[158];
    result_mem = num_digits_in_result * SIZES["digit"]

    # Variable 2: The main loop counter 'i' (from 2 to 100).
    # Max value is 100. 'cent' (0-99) is too small. 'char' (0-999) is required.
    # C declaration: char i;
    i_mem = SIZES["char"]

    # Variable 3: 'result_size', tracking the number of digits in the array.
    # Max value is 158. 'cent' (0-99) is too small. 'char' (0-999) is required.
    # C declaration: char result_size;
    result_size_mem = SIZES["char"]
    
    # Variable 4: Inner loop counter 'j'.
    # Max value is 157 (result_size - 1). 'cent' (0-99) is too small. 'char' (0-999) is required.
    # C declaration: char j;
    j_mem = SIZES["char"]

    # Variable 5: Temporary 'product' for multiplication step.
    # product = result[j] * i + carry. Max value is approx. 9 * 100 + 90 = 990.
    # This fits within a 'char' (0-999).
    # C declaration: char product;
    product_mem = SIZES["char"]

    # Variable 6: The 'carry' digit for multiplication.
    # carry = product / 10. Max value is 990 / 10 = 99.
    # This fits within a 'cent' (0-99).
    # C declaration: cent carry;
    carry_mem = SIZES["cent"]
    
    # The total memory 'z' is the sum of the memory for all variables.
    z = result_mem + i_mem + result_size_mem + j_mem + product_mem + carry_mem
    
    # Output the breakdown of the calculation for 'z'.
    print("--- Calculating z (Minimum Memory Size) ---")
    print(f"The calculation of 100! requires an array and several variables. Their memory sizes are:")
    print(f"1. Result array ('digit result[{num_digits_in_result}]'): {result_mem}D")
    print(f"2. Loop counter 'i' ('char'): {i_mem}D")
    print(f"3. Size tracker 'result_size' ('char'): {result_size_mem}D")
    print(f"4. Inner loop counter 'j' ('char'): {j_mem}D")
    print(f"5. Temp product ('char'): {product_mem}D")
    print(f"6. Carry variable ('cent'): {carry_mem}D")
    
    # Output the final equation for z as requested.
    print(f"\nFinal equation for z: {result_mem} + {i_mem} + {result_size_mem} + {j_mem} + {product_mem} + {carry_mem} = {z}")
    
    # Output the final answer in the format z:y
    print("\n--- Final Answer ---")
    print(f"{z}:{y}")

# Execute the solution
solve_wuxing_factorial()