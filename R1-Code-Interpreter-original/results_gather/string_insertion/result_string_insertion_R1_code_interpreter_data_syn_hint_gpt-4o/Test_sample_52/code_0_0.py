def apply_t_operation(s):
    i = 0
    result = []
    while i < len(s):
        if s[i:i+4] == "ABCD":
            result.append("ABCD")
            result.append("A")
            i += 4
        elif s[i:i+4] == "BCDE":
            result.append("BCDE")
            result.append("B")
            i += 4
        elif s[i:i+4] == "CDEA":
            result.append("CDEA")
            result.append("C")
            i += 4
        elif s[i:i+4] == "DEAB":
            result.append("DEAB")
            result.append("D")
            i += 4
        elif s[i:i+4] == "EABC":
            result.append("EABC")
            result.append("E")
            i += 4
        else:
            result.append(s[i])
            i += 1
    return ''.join(result)

s = "BBCDEBDEABADEABBA"
result = apply_t_operation(s)
print(result)