def operation1(s):
    return (s[:-1], True) if s.endswith('ca') else (s, False)

def operation2(s):
    if 'bca' in s:
        idx = s.find('bca')
        return (s[:idx] + s[idx+3:], True)
    return (s, False)

def operation3(s):
    return (s[:-2] + 'ab', True) if s.endswith('ba') else (s, False)

def operation4(s):
    return ('aa' + s[2:-1], True) if s.startswith('cb') else (s, False)

def operation5(s):
    if len(s) > 1 and 'ca' in s[1:]:
        idx = s[1:].find('ca') + 1
        return (s[:idx] + s[idx+2:], True)
    return (s, False)

def operation6(s):
    return (s[1:], True) if s.startswith('aa') else (s, False)

def process_string(s):
    print(f"Starting with: {s}")
    iteration = 1
    
    while True:
        original = s
        for op_num, op in enumerate([operation1, operation2, operation3, operation4, operation5, operation6], 1):
            s, changed = op(s)
            if changed:
                print(f"Iteration {iteration}, Operation {op_num}: {s}")
                break
        if s == original:
            break
        iteration += 1
    
    return s

final = process_string("bacbcbacc")
print(f"\nFinal result: {final}")