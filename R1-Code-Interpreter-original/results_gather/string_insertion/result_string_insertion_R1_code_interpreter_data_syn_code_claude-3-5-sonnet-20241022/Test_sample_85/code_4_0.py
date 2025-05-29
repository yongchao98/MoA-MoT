def apply_T(s):
    s = list(s)  # Convert to list for easier manipulation
    i = 0
    insertions = 0
    
    while i < len(s) - 3:  # We need at least 4 characters to check
        current = ''.join(s[i:i+4])  # Get current 4-char window
        
        if current == 'ABCD':
            s.insert(i + 4 + insertions, 'A')
            insertions += 1
        elif current == 'BCDE':
            s.insert(i + 4 + insertions, 'B')
            insertions += 1
        elif current == 'CDEA':
            s.insert(i + 4 + insertions, 'C')
            insertions += 1
        elif current == 'DEAB':
            s.insert(i + 4 + insertions, 'D')
            insertions += 1
        elif current == 'EABC':
            s.insert(i + 4 + insertions, 'E')
            insertions += 1
        i += 1
    
    result = ''.join(s)
    print(result)

input_string = "EABCBDEABADEAB"
apply_T(input_string)