def apply_T_operation(s):
    result = s
    i = 0
    while i <= len(result) - 4:
        substr = result[i:i+4]
        insert_char = None
        
        if substr == 'ABCD':
            insert_char = 'A'
        elif substr == 'BCDE':
            insert_char = 'B'
        elif substr == 'CDEA':
            insert_char = 'C'
        elif substr == 'DEAB':
            insert_char = 'D'
        elif substr == 'EABC':
            insert_char = 'E'
            
        if insert_char:
            result = result[:i+4] + insert_char + result[i+4:]
            i += 5  # Skip past the insertion
        else:
            i += 1
            
    print(result)

s = "BCDEEBCDECABCDDD"
apply_T_operation(s)