def apply_T_operation(s):
    # For debugging, let's print each match we find
    result = list(s)
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Find all matches and print them
    matches = []
    for i in range(len(s)):
        if i + 4 <= len(s):
            substr = s[i:i+4]
            if substr in patterns:
                matches.append((i, substr, patterns[substr]))
    
    # Print matches for verification
    for pos, pat, ins in matches:
        print(f"Found {pat} at position {pos}, will insert {ins}")
    
    # Apply insertions
    offset = 0
    for pos, pat, char in matches:
        insert_pos = pos + 4 + offset
        result.insert(insert_pos, char)
        offset += 1
    
    return ''.join(result)

s = "DEEABCCCDEAEBCDE"
result = apply_T_operation(s)
print("\nFinal result:", result)