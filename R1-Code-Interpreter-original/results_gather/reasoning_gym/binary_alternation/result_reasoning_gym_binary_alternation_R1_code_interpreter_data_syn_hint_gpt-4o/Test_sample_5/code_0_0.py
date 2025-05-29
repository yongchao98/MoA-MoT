def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to count mismatches for a given starting character
    def count_mismatches(start_char):
        mismatches = 0
        expected_char = start_char
        for char in s:
            if char != expected_char:
                mismatches += 1
            expected_char = '1' if expected_char == '0' else '0'
        return mismatches
    
    # Calculate swaps for both possible starting patterns
    if count_0 == count_1:
        swaps_start_with_0 = count_mismatches('0') // 2
        swaps_start_with_1 = count_mismatches('1') // 2
        return min(swaps_start_with_0, swaps_start_with_1)
    elif count_0 > count_1:
        return count_mismatches('0') // 2
    else:
        return count_mismatches('1') // 2

# Given binary string
binary_string = "10101010101111011010010100000"
# Calculate the minimum number of swaps
min_swaps = min_swaps_to_alternate(binary_string)
print(min_swaps)