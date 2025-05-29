def apply_T_operation(s):
    # Initialize result as list for easier manipulation
    result = list(s)
    offset = 0
    
    # Scan through the original string
    for i in range(len(s)-3):  # -3 to avoid index out of range
        # Check each pattern at current position
        substr = s[i:i+4]
        
        # Check for each pattern
        if substr == "ABCD":
            result.insert(i + 4 + offset, "A")
            offset += 1
        elif substr == "BCDE":
            result.insert(i + 4 + offset, "B")
            offset += 1
        elif substr == "CDEA":
            result.insert(i + 4 + offset, "C")
            offset += 1
        elif substr == "DEAB":
            result.insert(i + 4 + offset, "D")
            offset += 1
        elif substr == "EABC":
            result.insert(i + 4 + offset, "E")
            offset += 1
    
    final_string = ''.join(result)
    print(final_string)

# Test with the given string
s = "CBABCDADEABBABCD"
apply_T_operation(s)