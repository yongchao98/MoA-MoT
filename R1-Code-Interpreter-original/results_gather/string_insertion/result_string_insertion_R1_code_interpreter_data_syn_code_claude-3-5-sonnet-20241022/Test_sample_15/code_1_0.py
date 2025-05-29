def apply_T_operation(s):
    result = s
    i = 0
    while i < len(result):
        if i + 4 <= len(result):
            substr = result[i:i+4]
            if substr == "ABCD":
                result = result[:i+4] + "A" + result[i+4:]
                i += 5
            elif substr == "BCDE":
                result = result[:i+4] + "B" + result[i+4:]
                i += 5
            elif substr == "CDEA":
                result = result[:i+4] + "C" + result[i+4:]
                i += 5
            elif substr == "DEAB":
                result = result[:i+4] + "D" + result[i+4:]
                i += 5
            elif substr == "EABC":
                result = result[:i+4] + "E" + result[i+4:]
                i += 5
            else:
                i += 1
        else:
            i += 1
    return result

s = "ECDEAACDEAEDEAB"
print(apply_T_operation(s))