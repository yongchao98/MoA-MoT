def find_patterns(s):
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Store positions and characters to insert
    insertions = []
    
    # Scan string from left to right
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
    
    print(''.join(result))

# Test with given string
s = "CDEAECDEAAABCD"
find_patterns(s)