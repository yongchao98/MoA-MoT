def apply_T_operation(s):
    # Create a list to store positions and characters to insert
    insertions = []
    
    # Check each position for all possible patterns
    for i in range(len(s)):
        # Check ABCD pattern
        if i+3 < len(s) and s[i:i+4] == "ABCD":
            insertions.append((i+4, 'A'))
        # Check BCDE pattern
        if i+3 < len(s) and s[i:i+4] == "BCDE":
            insertions.append((i+4, 'B'))
        # Check CDEA pattern
        if i+3 < len(s) and s[i:i+4] == "CDEA":
            insertions.append((i+4, 'C'))
        # Check DEAB pattern
        if i+3 < len(s) and s[i:i+4] == "DEAB":
            insertions.append((i+4, 'D'))
        # Check EABC pattern
        if i+3 < len(s) and s[i:i+4] == "EABC":
            insertions.append((i+4, 'E'))
    
    # Sort insertions by position
    insertions.sort()
    
    # Apply insertions from left to right
    result = list(s)
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    print(''.join(result))
    print("Insertions:", insertions)

s = "ADEABBABCDAEABCBD"
apply_T_operation(s)