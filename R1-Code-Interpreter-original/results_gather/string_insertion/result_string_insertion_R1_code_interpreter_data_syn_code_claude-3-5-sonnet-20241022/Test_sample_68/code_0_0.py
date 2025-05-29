def find_all_matches(s):
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Store all matches with their positions
    matches = []
    for pattern, insert_char in patterns.items():
        pos = 0
        while True:
            pos = s.find(pattern, pos)
            if pos == -1:
                break
            matches.append((pos, pattern, insert_char))
            pos += 1
    
    # Sort matches by position
    matches.sort()
    
    if not matches:
        return s
    
    # Take the leftmost match
    pos, pattern, insert_char = matches[0]
    
    # Insert the character after the pattern
    result = s[:pos + len(pattern)] + insert_char + s[pos + len(pattern):]
    print(result)

# Test with the given string
s = "BCDEBDEABBDEAB"
find_all_matches(s)