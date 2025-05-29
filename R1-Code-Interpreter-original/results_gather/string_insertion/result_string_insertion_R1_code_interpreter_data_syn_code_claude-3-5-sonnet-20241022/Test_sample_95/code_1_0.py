def apply_T_operation(s):
    # Initialize result as list for easier manipulation
    result = list(s)
    offset = 0  # Track how many insertions we've made
    
    # Check each position from left to right
    for i in range(len(s)-3):  # -3 because we need 4 characters for each pattern
        # Get the current 4-character substring
        substr = s[i:i+4]
        
        # Check each pattern and insert accordingly
        if substr == 'ABCD':
            result.insert(i + 4 + offset, 'A')
            offset += 1
        elif substr == 'BCDE':
            result.insert(i + 4 + offset, 'B')
            offset += 1
        elif substr == 'CDEA':
            result.insert(i + 4 + offset, 'C')
            offset += 1
        elif substr == 'DEAB':
            result.insert(i + 4 + offset, 'D')
            offset += 1
        elif substr == 'EABC':
            result.insert(i + 4 + offset, 'E')
            offset += 1
    
    print(''.join(result))

# Test with the given string
s = "AADEABEDEABECDEA"
apply_T_operation(s)