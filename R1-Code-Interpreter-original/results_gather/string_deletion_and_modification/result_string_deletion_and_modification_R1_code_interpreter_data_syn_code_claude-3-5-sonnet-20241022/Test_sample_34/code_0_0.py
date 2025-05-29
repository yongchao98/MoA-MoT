def operation1(s):
    # If prefix is 'cb', replace with 'aa' and delete last char
    if s.startswith('cb'):
        return 'aa' + s[2:-1]
    return s

def operation2(s):
    # If suffix is 'bb', delete last two chars
    if s.endswith('bb'):
        return s[:-2]
    return s

def operation3(s):
    # If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        idx = s.index('acb')
        return s[:idx] + 'bca' + s[idx+3:]
    return s

def operation4(s):
    # If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa'
    return s

def operation5(s):
    # If starts with 'ac', replace first two chars with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:]
    return s

def operation6(s):
    # If contains 'ca' (not at start), remove first occurrence after first char
    if len(s) > 1:
        idx = s[1:].find('ca')
        if idx != -1:
            idx += 1  # Adjust index to account for the slice
            return s[:idx] + s[idx+2:]
    return s

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    
    print(f"Starting with: {s}")
    while True:
        original = s
        for op_num, op in enumerate(operations, 1):
            new_s = op(s)
            if new_s != s:
                print(f"Operation {op_num}: {s} -> {new_s}")
                s = new_s
                break
        if original == s:
            break
    return s

initial = "aacccabbbbcaacaacbcc"
final = process_string(initial)
print(f"\nFinal result: {final}")