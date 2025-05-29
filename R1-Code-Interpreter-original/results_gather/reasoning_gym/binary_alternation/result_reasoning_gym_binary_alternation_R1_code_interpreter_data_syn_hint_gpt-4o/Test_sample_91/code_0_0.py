def min_swaps_to_alternate(s):
    pattern1 = '01010101010'  # Starting with '0'
    pattern2 = '10101010101'  # Starting with '1'
    
    mismatches1 = sum(1 for i in range(len(s)) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(len(s)) if s[i] != pattern2[i])
    
    # Each swap can fix two mismatches
    swaps1 = mismatches1 // 2
    swaps2 = mismatches2 // 2
    
    return min(swaps1, swaps2)

s = "01001000111"
print(min_swaps_to_alternate(s))