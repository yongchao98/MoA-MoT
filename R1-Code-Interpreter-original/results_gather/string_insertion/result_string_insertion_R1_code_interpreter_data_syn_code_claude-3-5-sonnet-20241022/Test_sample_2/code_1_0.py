def find_all_patterns(s):
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Print each step for verification
    result = list(s)
    print("Original:", ''.join(result))
    
    for i in range(len(s)-3):
        substr = s[i:i+4]
        if substr in patterns:
            print(f"Found {substr} at position {i}")
            insert_pos = i + 4
            # Adjust the string
            result.insert(insert_pos, patterns[substr])
            s = ''.join(result)  # Update s for next iteration
            print("After insertion:", s)
    
    return ''.join(result)

s = "DACDEAADEABEEABCEB"
result = find_all_patterns(s)
print("\nFinal result:", result)