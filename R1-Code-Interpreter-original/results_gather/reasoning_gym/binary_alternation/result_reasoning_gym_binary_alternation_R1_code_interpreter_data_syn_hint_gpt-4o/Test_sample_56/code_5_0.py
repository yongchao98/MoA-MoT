def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to calculate swaps needed for a given pattern
    def calculate_swaps(start_char):
        swaps = 0
        for i in range(n):
            expected_char = start_char if i % 2 == 0 else ('1' if start_char == '0' else '0')
            if s[i] != expected_char:
                swaps += 1
        return swaps // 2
    
    # Calculate swaps for both patterns
    if count_0 == count_1:
        # If counts are equal, we can use either pattern
        swaps1 = calculate_swaps('0')
        swaps2 = calculate_swaps('1')
        return min(swaps1, swaps2)
    elif count_0 > count_1:
        # If there are more '0's, we must start with '0'
        return calculate_swaps('0')
    else:
        # If there are more '1's, we must start with '1'
        return calculate_swaps('1')

# Given binary string
binary_string = "010001001011111110100101"
result = min_swaps_to_alternate(binary_string)
print(result)