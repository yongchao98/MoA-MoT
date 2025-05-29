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
    
    # Scan the string
    for i in range(len(s)-3):  # -3 because patterns are 4 chars long
        substring = s[i:i+4]
        if substring in patterns:
            insertions.append((i+4, patterns[substring]))
    
    # Apply insertions from left to right
    # Sort by position to ensure left-to-right processing
    insertions.sort(key=lambda x: x[0])
    
    # Apply each insertion
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    return ''.join(result)

# Test the string
s = "EDDEABDABCDBBCDEC"
result = apply_T_operation(s)
print(result)