def apply_T_operation_with_steps(s):
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    result = list(s)
    changes = []
    
    # Find all matches and their positions
    for i in range(len(s)-3):
        substr = s[i:i+4]
        if substr in patterns:
            changes.append((i, substr, patterns[substr]))
    
    # Apply changes and show steps
    offset = 0
    for pos, pattern, insert_char in changes:
        result.insert(pos + 4 + offset, insert_char)
        offset += 1
        print(f"Found {pattern} at position {pos}, inserting {insert_char}")
        print(f"Current string: {''.join(result)}")
    
    print("\nFinal result:", ''.join(result))

s = "ACEABCBEABCBBCDEE"
apply_T_operation_with_steps(s)