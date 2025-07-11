import math

def solve_wuxing_factorial():
    """
    Calculates the memory cost (z) and first three digits (y)
    for computing 100! on the Wuxing virtual machine.
    """

    # --- Step 1: Define Wuxing data type sizes in decimal digits (D) ---
    CENT_SIZE = 2
    INT_SIZE = 6

    # --- Step 2: Calculate the optimized memory size (z) ---
    # The result of 100! has 157 decimal digits.
    # To store this, we use an array of 'cent' (0-99) to represent the number in base 100.
    # Number of array elements needed = ceil(157 / 2) = 79. We use 80 for safety.
    result_array_mem = 80 * CENT_SIZE

    # The main loop counter 'i' must go up to 100. A 'cent' (max 99) is too small.
    # Therefore, we must use an 'int'.
    loop_counter_i_mem = INT_SIZE

    # The inner loop counter 'j' and the array 'size' tracker will not exceed 80.
    # A 'cent' is sufficient for both.
    inner_loop_j_mem = CENT_SIZE
    size_tracker_mem = CENT_SIZE

    # The intermediate product can be up to `99 * 100 + carry` (~9999).
    # This requires an 'unsigned int' (max 999,999) for storage.
    product_mem = INT_SIZE
    carry_mem = INT_SIZE
    
    # The final equation for z is the sum of memory for all variables.
    z = result_array_mem + loop_counter_i_mem + inner_loop_j_mem + size_tracker_mem + product_mem + carry_mem
    
    print("Optimized Memory Calculation (z):")
    # This print statement fulfills the requirement: "output each number in the final equation"
    print(f"z = {result_array_mem} (result array) + {loop_counter_i_mem} (counter i) + {inner_loop_j_mem} (counter j) + {size_tracker_mem} (size tracker) + {product_mem} (product) + {carry_mem} (carry)")
    print(f"Total z = {z} D\n")


    # --- Step 3: Calculate the first three digits (y) of 100! ---
    factorial_100 = math.factorial(100)
    y = str(factorial_100)[:3]
    
    print("First 3 Digits of 100! (y):")
    print(f"100! starts with {y}...")
    print(f"y = {y}\n")
    
    # --- Step 4: Print the final answer in the specified z:y format ---
    print("Final Answer:")
    print(f"{z}:{y}")

solve_wuxing_factorial()