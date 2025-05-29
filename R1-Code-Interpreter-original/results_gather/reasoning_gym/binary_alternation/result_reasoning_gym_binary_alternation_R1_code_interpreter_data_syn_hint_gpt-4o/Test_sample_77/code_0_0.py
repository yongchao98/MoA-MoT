def min_swaps_to_alternate(s):
    n = len(s)
    
    # Pattern 1: "010101..."
    count1 = 0
    for i in range(n):
        expected_char = '0' if i % 2 == 0 else '1'
        if s[i] != expected_char:
            count1 += 1
    
    # Pattern 2: "101010..."
    count2 = 0
    for i in range(n):
        expected_char = '1' if i % 2 == 0 else '0'
        if s[i] != expected_char:
            count2 += 1
    
    # To make the string alternating, we need to swap half of the mismatches
    # because each swap fixes two positions.
    if count1 % 2 == 0 and count2 % 2 == 0:
        return min(count1 // 2, count2 // 2)
    elif count1 % 2 == 0:
        return count1 // 2
    elif count2 % 2 == 0:
        return count2 // 2
    else:
        return -1

# Given binary string
binary_string = "011110000001110110"
result = min_swaps_to_alternate(binary_string)
print(result)