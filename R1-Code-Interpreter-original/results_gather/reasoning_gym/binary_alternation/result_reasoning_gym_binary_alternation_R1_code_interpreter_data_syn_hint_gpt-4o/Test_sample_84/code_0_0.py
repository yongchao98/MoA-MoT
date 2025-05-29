def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to calculate swaps needed to match a given pattern
    def calculate_swaps(pattern):
        swaps = 0
        for i in range(n):
            if s[i] != pattern[i % 2]:
                swaps += 1
        return swaps // 2
    
    # Two possible patterns
    pattern1 = '01'  # Starting with '0'
    pattern2 = '10'  # Starting with '1'
    
    # Calculate swaps for both patterns
    if count_0 > count_1:
        # More '0's, must start with '0'
        return calculate_swaps(pattern1)
    elif count_1 > count_0:
        # More '1's, must start with '1'
        return calculate_swaps(pattern2)
    else:
        # Equal number of '0's and '1's, choose the minimum swaps
        return min(calculate_swaps(pattern1), calculate_swaps(pattern2))

# Given binary string
binary_string = "10011010100011111001110010100"
# Calculate the minimum swaps
min_swaps = min_swaps_to_alternate(binary_string)
print(min_swaps)