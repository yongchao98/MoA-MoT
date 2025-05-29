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
    
    # Scan the string from left to right
    for i in range(len(s)-3):  # -3 because patterns are 4 chars long
        substring = s[i:i+4]
        if substring in patterns:
            # Store the position where we'll insert and what character
            insertions.append((i+4, patterns[substring]))
    
    # Apply insertions from left to right
    # Need to adjust positions based on previous insertions
    offset = 0
    for pos, char in sorted(insertions):
        result.insert(pos + offset, char)
        offset += 1
    
    # Convert back to string
    result = ''.join(result)
    print(result)

# Test with the given string
s = "BCCDEAAEABCCCDEAD"
apply_T_operation(s)