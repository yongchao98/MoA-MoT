def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to calculate swaps needed for a given pattern
    def calculate_swaps(pattern):
        swaps = 0
        for i in range(n):
            if s[i] != pattern[i % 2]:
                swaps += 1
        return swaps // 2
    
    # Two possible patterns
    pattern1 = '01' * (n // 2) + '0' * (n % 2)
    pattern2 = '10' * (n // 2) + '1' * (n % 2)
    
    # Calculate swaps for both patterns
    swaps1 = calculate_swaps(pattern1)
    swaps2 = calculate_swaps(pattern2)
    
    # If the string length is even, both patterns are valid
    if n % 2 == 0:
        return min(swaps1, swaps2)
    else:
        # If the string length is odd, only one pattern is valid
        if count_0 > count_1:
            return swaps1
        else:
            return swaps2

# Given binary string
binary_string = "110010100010110001010101101001"
result = min_swaps_to_alternate(binary_string)
print(result)