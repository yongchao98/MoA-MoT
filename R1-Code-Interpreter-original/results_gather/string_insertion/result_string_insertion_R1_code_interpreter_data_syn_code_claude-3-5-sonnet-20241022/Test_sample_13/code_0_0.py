def operation_T(s):
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
    insertions = []  # Store (position, char) pairs for insertion
    
    # Scan the string from left to right
    for i in range(len(s)-3):
        substring = s[i:i+4]
        if substring in patterns:
            # Store position and character to insert
            insertions.append((i+4, patterns[substring]))
    
    # Apply insertions from left to right
    # Need to adjust positions based on previous insertions
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    final_string = ''.join(result)
    print(final_string)

# Test with the given string
s = "DEABDDEABADEABBA"
operation_T(s)