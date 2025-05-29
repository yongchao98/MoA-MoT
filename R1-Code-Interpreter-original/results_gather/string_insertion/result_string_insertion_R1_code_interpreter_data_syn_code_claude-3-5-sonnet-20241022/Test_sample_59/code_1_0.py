def apply_T_operation(s):
    # Convert string to list for easier manipulation
    result = list(s)
    
    # Store all positions where we need to insert characters
    insertions = []
    
    # Check each position in the string
    i = 0
    while i <= len(s) - 4:
        substr = s[i:i+4]
        insert_char = None
        
        if substr == "ABCD":
            insert_char = "A"
        elif substr == "BCDE":
            insert_char = "B"
        elif substr == "CDEA":
            insert_char = "C"
        elif substr == "DEAB":
            insert_char = "D"
        elif substr == "EABC":
            insert_char = "E"
            
        if insert_char:
            insertions.append((i+4, insert_char))
            i += 4  # Skip to position after the pattern
        else:
            i += 1
    
    # Apply insertions from left to right
    offset = 0
    for pos, char in insertions:
        result.insert(pos + offset, char)
        offset += 1
    
    final_result = ''.join(result)
    print(final_result)

# Test with the given string
s = "ADEABEDEABBEABCBE"
apply_T_operation(s)