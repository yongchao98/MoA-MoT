def apply_T_operation(s):
    # Define patterns and their corresponding insertions
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Keep track of insertions (position, character)
    insertions = []
    
    # Scan string from left to right
    for i in range(len(s)-3):  # -3 because patterns are 4 chars long
        substring = s[i:i+4]
        if substring in patterns:
            # When pattern found, record position and character to insert
            insertions.append((i+4, patterns[substring]))
    
    # Convert string to list for easier insertion
    result = list(s)
    
    # Apply insertions from left to right
    # Need to adjust positions based on previous insertions
    offset = 0
    for pos, char in sorted(insertions):
        result.insert(pos + offset, char)
        offset += 1
    
    # Convert back to string
    result = ''.join(result)
    print(result)

# Test with given string
s = "EDBCDEEBCDECABCDCB"
apply_T_operation(s)