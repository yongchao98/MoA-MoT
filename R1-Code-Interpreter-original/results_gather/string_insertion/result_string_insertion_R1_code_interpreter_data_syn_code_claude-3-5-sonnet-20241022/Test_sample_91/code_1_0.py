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
    i = 0
    while i <= len(s)-4:  # Check all possible 4-char windows
        substring = s[i:i+4]
        if substring in patterns:
            insertions.append((i+4, patterns[substring]))
            i += 4  # Move past the current pattern
        else:
            i += 1
    
    # Apply insertions from left to right
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    final_result = ''.join(result)
    print(f"Original: {s}")
    print(f"Final: {final_result}")
    print("Patterns found:")
    i = 0
    while i <= len(s)-4:
        substring = s[i:i+4]
        if substring in patterns:
            print(f"Position {i}: {substring} -> insert {patterns[substring]}")
            i += 4
        else:
            i += 1

# Test with the given string
s = "AECDEACCDEAACDEAAB"
apply_T_operation(s)