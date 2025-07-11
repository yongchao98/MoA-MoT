import math

def solve_wuxing_factorial():
    """
    This function calculates the memory footprint and first three digits of 100!
    based on the constraints of the Wuxing virtual machine.
    """

    # Part 1: Calculate the smallest memory size (z) in Decimal Digits (D)
    # -------------------------------------------------------------------
    # Strategy: Use an array of `char` (3D, 0-999) to store 100! in base 1000.
    # This is more memory-efficient than a simple array of `digit`s (base 10).

    # 100! has 158 decimal digits (math.log10(math.factorial(100)) is approx 157.97).
    # Number of base-1000 chunks needed = ceil(158 / 3) = 53 chunks.
    num_chunks = 53

    # Define Wuxing data type sizes in Decimal Digits (D)
    D_CHAR = 3
    D_CENT = 2
    D_UINT = 6

    # Determine variable types and calculate total memory size 'z'
    # 1. `result` array: Stores the large number. Needs 53 elements.
    size_result = num_chunks * D_CHAR  # char result[53];

    # 2. Loop counter `i` (from 2 to 100): `cent` (0-99) is too small. `char` (0-999) is needed.
    size_i = D_CHAR  # char i;

    # 3. Inner loop counter `j` (from 0 to 52): `cent` (0-99) is sufficient.
    size_j = D_CENT  # cent j;

    # 4. `num_chunks` variable (tracks active chunks, max 53): `cent` (0-99) is sufficient.
    size_num_chunks = D_CENT # cent num_chunks;

    # 5. `carry` in multiplication: Max value is (999 * 100 + 99) / 1000 = 99. `cent` is sufficient.
    size_carry = D_CENT  # cent carry;

    # 6. `product` variable: Max value is 999 * 100 + 99 = 99999. Needs 6 digits.
    #    `unsigned int` (0-999,999) is the smallest type that can hold this.
    size_product = D_UINT  # unsigned int product;

    # Total memory size 'z'
    z = size_result + size_i + size_j + size_num_chunks + size_carry + size_product

    # Part 2: Calculate the first 3 digits of 100! (y)
    # ----------------------------------------------------
    # We simulate the base-1000 calculation.
    BASE = 1000
    
    # Initialize result to 1
    res_array = [1]
    
    # Loop from i = 2 to 100
    for i in range(2, 101):
        carry = 0
        for j in range(len(res_array)):
            product = res_array[j] * i + carry
            res_array[j] = product % BASE
            carry = product // BASE
        
        while carry > 0:
            res_array.append(carry % BASE)
            carry //= BASE
            
    # The result is stored with the least significant chunk first. We reverse it.
    res_array.reverse()
    
    # The first chunk (res_array[0]) contains the most significant digits.
    most_significant_chunk = res_array[0]
    
    # The next chunk is needed in case the first one has fewer than 3 digits.
    next_significant_chunk = res_array[1]
    
    # Combine the first two chunks and take the first 3 digits.
    # Format the second chunk to be 3 digits long with leading zeros if necessary.
    full_number_str = str(most_significant_chunk) + f"{next_significant_chunk:03d}"
    y = full_number_str[:3]

    # Part 3: Print the final answer in the required format z:y
    # ------------------------------------------------------------
    # The prompt requires outputting each number in the final equation.
    # We will print z and y individually before the final formatted output.
    print(f"Calculated smallest memory size (z): {z}")
    print(f"Calculated first 3 digits of 100! (y): {y}")
    print(f"{z}:{y}")


solve_wuxing_factorial()