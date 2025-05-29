def operation_T(s):
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
    insertions = []  # Store insertions to be made
    
    # Scan the string
    while i <= len(result) - 4:
        substring = ''.join(result[i:i+4])
        if substring in patterns:
            # Store the position and character to insert
            insertions.append((i+4, patterns[substring]))
        i += 1
    
    # Apply insertions from left to right
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    return ''.join(result)

# Test with the given string
s = "ADABCDCCDEADABCD"
result = operation_T(s)
print(result)