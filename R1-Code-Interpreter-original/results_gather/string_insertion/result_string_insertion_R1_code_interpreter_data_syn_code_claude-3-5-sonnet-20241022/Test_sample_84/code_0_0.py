def find_all_patterns(s):
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Store positions and characters to insert
    insertions = []
    
    # Check each position
    for i in range(len(s)-3):
        substr = s[i:i+4]
        if substr in patterns:
            insertions.append((i+4, patterns[substr]))
    
    # Sort insertions by position
    insertions.sort()
    
    # Create new string with insertions
    result = list(s)
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    return ''.join(result)

s = "BBCDECBCDEECDEA"
result = find_all_patterns(s)
print(result)