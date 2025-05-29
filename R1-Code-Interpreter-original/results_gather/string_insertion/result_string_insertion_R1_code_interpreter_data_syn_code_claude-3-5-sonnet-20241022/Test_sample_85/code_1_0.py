def apply_T(s):
    # Original string length
    n = len(s)
    result = []
    i = 0
    
    # Process the string character by character
    while i < n-3:  # Need at least 4 characters ahead
        result.append(s[i])
        substr = s[i:i+4]
        
        # Check each pattern once
        if substr == 'ABCD':
            result.append(s[i+1])
            result.append(s[i+2])
            result.append(s[i+3])
            result.append('A')
            i += 4
        elif substr == 'BCDE':
            result.append(s[i+1])
            result.append(s[i+2])
            result.append(s[i+3])
            result.append('B')
            i += 4
        elif substr == 'CDEA':
            result.append(s[i+1])
            result.append(s[i+2])
            result.append(s[i+3])
            result.append('C')
            i += 4
        elif substr == 'DEAB':
            result.append(s[i+1])
            result.append(s[i+2])
            result.append(s[i+3])
            result.append('D')
            i += 4
        elif substr == 'EABC':
            result.append(s[i+1])
            result.append(s[i+2])
            result.append(s[i+3])
            result.append('E')
            i += 4
        else:
            i += 1
    
    # Append remaining characters
    result.extend(s[i:])
    print(''.join(result))

input_string = "EABCBDEABADEAB"
apply_T(input_string)