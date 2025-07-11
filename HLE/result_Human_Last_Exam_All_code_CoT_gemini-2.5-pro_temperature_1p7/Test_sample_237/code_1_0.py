def solve_factorial_on_wuxing():
    """
    This script calculates 100! using a BigInt algorithm suited for the
    hypothetical XVM architecture. It determines the minimum memory usage 'z'
    and the first three digits 'y' of the result.
    """

    # --- Step 1: BigInt Calculation for 100! ---
    # We represent the number as a list of integers, where each integer is a
    # chunk of the number in base 1000. This simulates an array of 'char' (3D, 0-999) on XVM.
    # The list is stored in reverse order (least significant chunk first).
    # For example, the number 9,876,543 would be stored as [543, 876, 9].
    result_b1000 = [1]
    
    # We track the maximum length of the list to determine the required array size.
    max_len_result = 1

    # Loop from i = 2 to 100 to compute the factorial.
    for i in range(2, 101):
        # The 'carry' can reach up to 100 during multiplication, so a 'char' (3D) is the
        # smallest sufficient type. A 'cent' (2D, 0-99) would be too small.
        carry = 0

        # Iterate through each chunk of our number.
        for j in range(len(result_b1000)):
            # The 'product' must hold 'chunk * i + carry'. The max value is roughly
            # 999 * 100 + 100 = 100,000. An 'unsigned int' (6D, 0-999,999) is needed.
            product = result_b1000[j] * i + carry
            result_b1000[j] = product % 1000
            carry = product // 1000
        
        # If there's a remaining carry, append it as new chunks.
        while carry > 0:
            result_b1000.append(carry % 1000)
            carry //= 1000
        
        # Update the tracker for the maximum list length.
        if len(result_b1000) > max_len_result:
            max_len_result = len(result_b1000)

    # --- Step 2: Calculate 'z' (Smallest Memory Size in D) ---
    # 'result[]': The array to hold the number. Its size is 'max_len_result'.
    # Each element is a 'char' (3D). max_len_result will be 53 for 100!.
    # z_result = 53 * 3D = 159D.
    # 'size': A variable to track the active size of the result array. max_len_result is ~53.
    # A 'cent' (2D) is sufficient.
    # 'i': Loop counter up to 100. A 'char' (3D) is the smallest sufficient type.
    # 'j': Inner loop counter up to 'size'. A 'cent' (2D) is sufficient.
    # 'carry': As analyzed, needs a 'char' (3D).
    # 'product': As analyzed, needs an 'unsigned int' (6D).
    # Total z = z_result + z_size + z_i + z_j + z_carry + z_product
    # z = (53 * 3) + 2 + 3 + 2 + 3 + 6 = 159 + 16 = 175.

    # --- Step 3: Calculate 'y' (First 3 Digits) and Print Full Result ---
    # The result is stored in reverse. The most significant chunk is the last element.
    most_significant_chunk = result_b1000[-1]
    second_most_significant_chunk = result_b1000[-2]
    
    # Reconstruct the start of the number to find the first 3 digits.
    # e.g., if chunks are [..., 332, 93], the number is 93,332,...
    # The string representation starts with "93" followed by "332" (padded).
    y = (str(most_significant_chunk) + f'{second_most_significant_chunk:03d}')[:3]

    # As per the instruction "output each number in the final equation!", we print the full 100!.
    # Start with the most significant chunk (no padding needed).
    final_number_str = str(result_b1000[-1])
    
    # Append the rest of the chunks, padding with leading zeros to 3 digits.
    for k in range(len(result_b1000) - 2, -1, -1):
        final_number_str += f'{result_b1000[k]:03d}'
    
    print(final_number_str)

solve_factorial_on_wuxing()