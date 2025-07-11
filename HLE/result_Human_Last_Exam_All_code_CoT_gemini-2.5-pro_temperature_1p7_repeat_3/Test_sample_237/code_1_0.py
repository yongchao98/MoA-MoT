import math

def solve_wuxing_factorial():
    """
    Solves the Wuxing factorial problem by calculating memory usage (z)
    and the first three digits of the result (y).
    """

    # Step 1: Analyze the properties of 100!
    factorial_100 = math.factorial(100)
    num_digits = len(str(factorial_100))

    # Step 2: Determine the optimal storage method.
    # We will use an array of 'char' types (3D, base 1000) as it's the most memory-efficient.
    # Number of char elements needed for the result array.
    num_chars_in_array = math.ceil(num_digits / 3)

    # Step 3: Calculate the memory size (z) for all variables in the C program.
    # Data type sizes in Decimal Digits (D)
    DIGIT_SIZE = 1
    CENT_SIZE = 2
    CHAR_SIZE = 3
    INT_SIZE = 6

    # Memory for the array holding the result (res)
    size_res_array = num_chars_in_array * CHAR_SIZE

    # Memory for helper variables
    # res_size (stores up to 53) -> cent
    size_res_size = CENT_SIZE
    # Loop counter i (up to 100) -> cent
    size_i = CENT_SIZE
    # Loop counter j (up to 52) -> cent
    size_j = CENT_SIZE
    # Carry variable (can hold values up to ~100,000) -> int
    size_carry = INT_SIZE
    # Printing variables (d1, d2, d3 to hold single digits) -> digit * 3
    size_print_vars = 3 * DIGIT_SIZE

    # Sum of all memory allocations gives z.
    z = size_res_array + size_res_size + size_i + size_j + size_carry + size_print_vars

    # Step 4: Determine the first 3 digits (y) of 100!
    y = int(str(factorial_100)[:3])

    # Step 5: Print the detailed calculation and the final result as requested.
    print("--- Wuxing 100! Calculation Analysis ---")
    print(f"Value of 100! has {num_digits} digits.")
    print(f"Optimal storage is an array of 'char' (3D), requiring {num_chars_in_array} elements.\n")

    print("--- Minimum Memory Calculation (z) ---")
    print("Variable    | Wuxing Type | Memory (D) | Justification")
    print("------------|-------------|------------|-----------------------------------------")
    print(f"res[{num_chars_in_array}]   | char[]      | {size_res_array:<10} | To store {num_digits}-digit result in 3D chunks")
    print(f"res_size    | cent        | {size_res_size:<10} | To track array size (up to {num_chars_in_array})")
    print(f"i           | cent        | {size_i:<10} | Loop counter for factorial (2 to 100)")
    print(f"j           | cent        | {size_j:<10} | Inner loop counter (up to {num_chars_in_array-1})")
    print(f"carry       | int         | {size_carry:<10} | To handle product/carry during multiplication")
    print(f"d1, d2, d3  | digit[3]    | {size_print_vars:<10} | To print first 3 digits individually")
    print("------------|-------------|------------|-----------------------------------------")
    
    # "output each number in the final equation!"
    print(f"\nFinal Equation for z: {z} = {size_res_array} + {size_res_size} + {size_i} + {size_j} + {size_carry} + {size_print_vars}")

    print(f"\n--- Final Result ---")
    print(f"z (Smallest Memory Size): {z} D")
    print(f"y (First 3 digits of 100!): {y}")
    print("\nFormatted Answer:")
    print(f"{z}:{y}")

solve_wuxing_factorial()