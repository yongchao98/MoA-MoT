def min_swaps_to_alternate(s):
    n = len(s)
    pattern1 = '01' * (n // 2) + '0' * (n % 2)
    pattern2 = '10' * (n // 2) + '1' * (n % 2)
    
    mismatches1 = sum(1 for i in range(n) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(n) if s[i] != pattern2[i])
    
    if n % 2 == 0:
        # If length is even, both patterns are possible
        return min(mismatches1 // 2, mismatches2 // 2)
    else:
        # If length is odd, only one pattern is possible
        if mismatches1 % 2 == 0:
            return mismatches1 // 2
        elif mismatches2 % 2 == 0:
            return mismatches2 // 2
        else:
            return -1

s = "1011000000110101001011110"
print(min_swaps_to_alternate(s))