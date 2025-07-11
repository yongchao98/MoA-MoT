import math

def calculate_factorial_on_xvm():
    """
    Simulates the calculation of 100! on the XVM to find the minimum
    memory usage (z) and the first 3 digits of the result (y).
    """

    # 1. Define memory costs for XVM data types in decimal digits (D)
    memory_costs = {
        'cent': 2,
        'char': 3,
        'int': 6,
    }

    # 2. Define variables and calculate memory cost 'z'
    # We store the 158-digit result in an array of 2-digit chunks ('cent').
    # Array size needed = 158 / 2 = 79. We use 80 for safety.
    ARRAY_SIZE = 80
    var_result_cost = ARRAY_SIZE * memory_costs['cent']

    # Loop counter i (2 to 100): 'char' (0-999) is the smallest type.
    var_i_cost = memory_costs['char']

    # Loop counter j (0 to 79): 'cent' (0-99) is sufficient.
    var_j_cost = memory_costs['cent']

    # num_digits (tracks array size): Max is ~79, 'cent' is sufficient.
    var_num_digits_cost = memory_costs['cent']
    
    # Intermediate multiplication and carry can exceed 99*100, so we need 'int'.
    var_temp_cost = memory_costs['int']
    var_carry_cost = memory_costs['int']
    
    # Final 3-digit result 'y' fits in a 'char'.
    var_y_cost = memory_costs['char']

    # Total memory cost 'z'
    z = (var_result_cost + var_i_cost + var_j_cost +
         var_num_digits_cost + var_temp_cost + var_carry_cost + var_y_cost)

    # 3. Implement the large number factorial algorithm
    # Each element in the list represents a 2-digit chunk (0-99).
    BASE = 100
    
    # 'result' array holds the number, stored in reverse order (least significant chunk first).
    result = [0] * ARRAY_SIZE
    result[0] = 1
    num_digits = 1

    # Multiply by i from 2 to 100
    for i in range(2, 101):
        carry = 0
        # For each 2-digit chunk in our number
        for j in range(num_digits):
            # 'temp' stores the result of multiplying the chunk by i and adding the carry.
            # This is the core of the "by-hand" multiplication algorithm.
            temp = result[j] * i + carry
            result[j] = temp % BASE
            carry = temp // BASE
        
        # If there's a carry left over, add it as new chunks to the end of the array.
        while carry > 0:
            result[num_digits] = carry % BASE
            carry = carry // BASE
            num_digits += 1

    # 4. Extract the first 3 digits 'y'
    # The most significant chunks are at the end of the array.
    last_chunk = result[num_digits - 1]
    second_to_last_chunk = result[num_digits - 2]
    
    # The last chunk (e.g., 93) gives the first two digits.
    # The first digit of the second-to-last chunk (e.g., 3 from 32) gives the third digit.
    y = last_chunk * 10 + (second_to_last_chunk // 10)

    # 5. Print the final answer in the format z:y
    print(f"Minimal memory z = {z} D")
    print(f"First 3 digits of 100! = {y}")
    print(f"Final Answer: {z}:{y}")

calculate_factorial_on_xvm()