def min_swaps_to_alternate(s):
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(len(s))])
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(len(s))])
    
    mismatches1 = sum(1 for i in range(len(s)) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(len(s)) if s[i] != pattern2[i])
    
    return min(mismatches1, mismatches2)

s = "11010010001"
print(min_swaps_to_alternate(s))