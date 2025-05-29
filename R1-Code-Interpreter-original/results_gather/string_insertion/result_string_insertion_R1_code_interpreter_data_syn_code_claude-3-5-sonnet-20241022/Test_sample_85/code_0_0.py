def apply_T(s):
    patterns = {
        'ABCD': 'A',
        'BCDE': 'B',
        'CDEA': 'C',
        'DEAB': 'D',
        'EABC': 'E'
    }
    
    # Convert to list for easier manipulation
    s = list(s)
    i = 0
    
    while i <= len(s) - 4:  # Need at least 4 characters to check
        substring = ''.join(s[i:i+4])
        found = False
        
        for pattern, insert_char in patterns.items():
            if substring == pattern:
                s.insert(i + 4, insert_char)
                found = True
                break
                
        i += 1 if not found else 0
    
    return ''.join(s)

input_string = "EABCBDEABADEAB"
result = apply_T(input_string)
print(result)