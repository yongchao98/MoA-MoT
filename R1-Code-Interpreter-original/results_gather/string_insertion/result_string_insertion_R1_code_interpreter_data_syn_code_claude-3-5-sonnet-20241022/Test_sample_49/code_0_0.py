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
    i = 0
    insertions_made = 0
    
    while i <= len(s) - 4:
        substring = s[i:i+4]
        if substring in patterns:
            # Insert the corresponding character at the correct position
            insert_pos = i + 4 + insertions_made
            result.insert(insert_pos, patterns[substring])
            insertions_made += 1
        i += 1
    
    return ''.join(result)

# Test with the given string
s = "AAABCDBDEABADEAB"
result = apply_T_operation(s)
print(result)