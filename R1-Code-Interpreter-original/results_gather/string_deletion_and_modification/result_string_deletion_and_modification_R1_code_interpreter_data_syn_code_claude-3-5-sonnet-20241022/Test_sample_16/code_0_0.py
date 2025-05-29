def operation1(s):
    # If suffix is 'bb', delete last two characters
    if s.endswith('bb'):
        return s[:-2], True
    return s, False

def operation2(s):
    # If contains 'abc', replace first occurrence with 'cab'
    if 'abc' in s:
        idx = s.index('abc')
        return s[:idx] + 'cab' + s[idx+3:], True
    return s, False

def operation3(s):
    # If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    return s, False

def operation4(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def operation5(s):
    # If contains 'ca' (not at start), remove first occurrence after first char
    if 'ca' in s[1:]:
        idx = s[1:].index('ca') + 1
        return s[:idx] + s[idx+2:], True
    return s, False

def operation6(s):
    # If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    return s, False

def process_string(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    changes_made = True
    steps = []
    
    while changes_made:
        changes_made = False
        for op in operations:
            new_s, changed = op(s)
            if changed:
                steps.append(f"{s} -> {new_s}")
                s = new_s
                changes_made = True
                break
    
    print("\n".join(steps))
    print("\nFinal string:", s)

# Test with the given string
process_string("bbabccbacbacbacbaaa")