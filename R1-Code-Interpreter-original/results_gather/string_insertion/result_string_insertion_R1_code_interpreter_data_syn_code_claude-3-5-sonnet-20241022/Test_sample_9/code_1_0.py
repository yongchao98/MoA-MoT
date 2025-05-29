def apply_T_operation(s):
    # Initialize result string
    result = list(s)
    offset = 0
    
    # Scan through original string
    for i in range(len(s)):
        # Check each pattern at current position
        if i+4 <= len(s):
            substr = s[i:i+4]
            insert_pos = i + 4 + offset
            
            if substr == "ABCD":
                result.insert(insert_pos, "A")
                offset += 1
            elif substr == "BCDE":
                result.insert(insert_pos, "B")
                offset += 1
            elif substr == "CDEA":
                result.insert(insert_pos, "C")
                offset += 1
            elif substr == "DEAB":
                result.insert(insert_pos, "D")
                offset += 1
            elif substr == "EABC":
                result.insert(insert_pos, "E")
                offset += 1
    
    return ''.join(result)

s = "EAEABCBDEABACDEADB"
print(apply_T_operation(s))