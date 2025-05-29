def find_all_patterns(s):
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Store (position, character) pairs for insertions
    insertions = []
    
    # Check each position
    for i in range(len(s)):
        # Get substring of length 4 starting at position i
        substr = s[i:i+4]
        if substr in patterns:
            # Record the position after the pattern and the character to insert
            insertions.append((i+4, patterns[substr]))
    
    # Sort insertions by position (though they should already be in order)
    insertions.sort(key=lambda x: x[0])
    
    # Create new string with insertions
    result = list(s)
    # Keep track of how many insertions we've made to adjust positions
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    return ''.join(result)

# Test with the given string
s = "BCDECCDEAACDEA"
result = find_all_patterns(s)
print(result)