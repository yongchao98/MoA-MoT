def apply_T(s):
    result = []
    i = 0
    
    while i < len(s):
        result.append(s[i])
        if i <= len(s) - 4:
            substr = s[i:i+4]
            if substr == 'ABCD':
                if i+4 < len(s):
                    result.append(s[i+1:i+4])
                result.append('A')
                i += 3
            elif substr == 'BCDE':
                if i+4 < len(s):
                    result.append(s[i+1:i+4])
                result.append('B')
                i += 3
            elif substr == 'CDEA':
                if i+4 < len(s):
                    result.append(s[i+1:i+4])
                result.append('C')
                i += 3
            elif substr == 'DEAB':
                if i+4 < len(s):
                    result.append(s[i+1:i+4])
                result.append('D')
                i += 3
            elif substr == 'EABC':
                if i+4 < len(s):
                    result.append(s[i+1:i+4])
                result.append('E')
                i += 3
        i += 1
    
    final = ''.join(result)
    print(final)

s = "EABCBDEABADEAB"
apply_T(s)