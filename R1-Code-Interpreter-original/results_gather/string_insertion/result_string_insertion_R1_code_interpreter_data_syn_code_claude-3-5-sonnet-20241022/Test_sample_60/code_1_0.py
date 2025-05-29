def apply_T_operation(s):
    # Define the patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Convert to list for easier manipulation
    result = list(s)
    insertions = []  # Store (position, char) pairs
    
    # Print each substring and pattern found for verification
    for i in range(len(s)):
        if i + 4 <= len(s):
            substr = s[i:i+4]
            if substr in patterns:
                print(f"Found pattern {substr} at position {i}")
                insertions.append((i+4, patterns[substr]))
    
    # Apply insertions from left to right
    offset = 0
    for pos, char in sorted(insertions):
        result.insert(pos + offset, char)
        offset += 1
    
    final_result = ''.join(result)
    print("Patterns found:", insertions)
    return final_result

# Test with the given string
s = "BCDEDDEABBEABCAE"
result = apply_T_operation(s)
print("Final result:", result)