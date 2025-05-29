def apply_T(s):
    # Convert string to list for easier manipulation
    result = list(s)
    offset = 0
    
    # Process string from left to right
    for i in range(len(s)-3):
        # Get current 4-char window
        window = s[i:i+4]
        # Check each pattern
        if window == 'ABCD':
            result.insert(i + 4 + offset, 'A')
            offset += 1
        elif window == 'BCDE':
            result.insert(i + 4 + offset, 'B')
            offset += 1
        elif window == 'CDEA':
            result.insert(i + 4 + offset, 'C')
            offset += 1
        elif window == 'DEAB':
            result.insert(i + 4 + offset, 'D')
            offset += 1
        elif window == 'EABC':
            result.insert(i + 4 + offset, 'E')
            offset += 1
    
    return ''.join(result)

s = "ABCDDCDEADBCDE"
print(apply_T(s))