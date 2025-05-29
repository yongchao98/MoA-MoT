def apply_T_operation(s):
    # Define the patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    result = list(s)
    offset = 0
    
    # Scan through original string positions
    for i in range(len(s)):
        if i + 4 <= len(s):
            substr = s[i:i+4]
            # Print each check for verification
            if substr in patterns:
                print(f"Found {substr} at position {i}, inserting {patterns[substr]}")
                result.insert(i + 4 + offset, patterns[substr])
                offset += 1
    
    final_result = ''.join(result)
    print(f"\nFinal result: {final_result}")

s = "EBCDEEEABCCBCDE"
apply_T_operation(s)