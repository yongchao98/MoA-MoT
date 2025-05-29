def apply_T_operation(s):
    # Store the original string length
    original_length = len(s)
    
    # Convert to list for easier manipulation
    result = list(s)
    insertions = []
    
    # Scan for all patterns from left to right
    i = 0
    while i <= len(s) - 4:
        substr = s[i:i+4]
        if substr == 'ABCD':
            insertions.append((i+4, 'A'))
        elif substr == 'BCDE':
            insertions.append((i+4, 'B'))
        elif substr == 'CDEA':
            insertions.append((i+4, 'C'))
        elif substr == 'DEAB':
            insertions.append((i+4, 'D'))
        elif substr == 'EABC':
            insertions.append((i+4, 'E'))
        i += 1
    
    # Sort insertions by position
    insertions.sort()
    
    # Apply insertions
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    final_string = ''.join(result)
    print(f"Original: {s}")
    print(f"Final: {final_string}")
    print(f"Insertions: {insertions}")

s = "BBEABCCABCDADEABB"
apply_T_operation(s)