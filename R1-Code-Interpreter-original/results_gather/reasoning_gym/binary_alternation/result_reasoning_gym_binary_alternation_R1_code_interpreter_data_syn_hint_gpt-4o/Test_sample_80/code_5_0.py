def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between counts of '0's and '1's is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate mismatches for pattern starting with '0' (010101...)
    mismatches_0 = sum(1 for i in range(n) if s[i] != str(i % 2))
    
    # Calculate mismatches for pattern starting with '1' (101010...)
    mismatches_1 = sum(1 for i in range(n) if s[i] != str((i + 1) % 2))
    
    # If the string length is even, both patterns are possible
    if n % 2 == 0:
        return min(mismatches_0, mismatches_1) // 2
    else:
        # If the string length is odd, only one pattern is possible
        if count_0 > count_1:
            return mismatches_0 // 2
        else:
            return mismatches_1 // 2

# Given binary string
binary_string = "00101010110010110011100001"
result = min_swaps_to_alternate(binary_string)
print(result)