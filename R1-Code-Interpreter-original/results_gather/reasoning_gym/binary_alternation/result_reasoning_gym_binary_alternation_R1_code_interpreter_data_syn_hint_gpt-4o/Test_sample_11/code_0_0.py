def min_swaps_to_alternate(s):
    pattern1 = '01' * (len(s) // 2) + '0' * (len(s) % 2)
    pattern2 = '10' * (len(s) // 2) + '1' * (len(s) % 2)
    
    mismatches1 = sum(1 for i in range(len(s)) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(len(s)) if s[i] != pattern2[i])
    
    # Each swap can fix two mismatches
    swaps1 = mismatches1 // 2
    swaps2 = mismatches2 // 2
    
    return min(swaps1, swaps2)

s = "00110110100"
print(min_swaps_to_alternate(s))