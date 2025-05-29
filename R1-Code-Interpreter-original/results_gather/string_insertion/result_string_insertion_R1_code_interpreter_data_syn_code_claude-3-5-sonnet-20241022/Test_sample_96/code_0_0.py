def apply_T(s):
    # Define the patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Convert string to list for easier manipulation
    result = list(s)
    # Keep track of how many characters we've inserted
    offset = 0
    
    # Check each position in original string
    for i in range(len(s)-3):  # -3 because patterns are 4 chars long
        # Get the 4-character substring starting at position i
        substr = s[i:i+4]
        # Check if this substring matches any pattern
        if substr in patterns:
            # Insert the corresponding character after the pattern
            insert_pos = i + 4 + offset
            result.insert(insert_pos, patterns[substr])
            offset += 1
    
    # Convert back to string and return
    return ''.join(result)

# Test with the given string
s = "ABCDDCDEADBCDE"
print(apply_T(s))