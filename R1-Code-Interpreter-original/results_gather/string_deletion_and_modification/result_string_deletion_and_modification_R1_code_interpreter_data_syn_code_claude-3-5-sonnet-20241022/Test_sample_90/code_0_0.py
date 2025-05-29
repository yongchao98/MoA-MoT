def operation1(s):  # length > 15, remove middle char
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def operation2(s):  # starts with 'bb', remove second char
    if s.startswith('bb'):
        return s[0] + s[2:], True
    return s, False

def operation3(s):  # suffix 'ac', replace with 'cb'
    if s.endswith('ac'):
        return s[:-2] + 'cb', True
    return s, False

def operation4(s):  # contains 'ca' not at start, remove first such occurrence
    if len(s) > 1 and 'ca' in s[1:]:
        idx = s.find('ca', 1)
        return s[:idx] + s[idx+2:], True
    return s, False

def operation5(s):  # contains 'abc', replace with 'cab'
    if 'abc' in s:
        idx = s.find('abc')
        return s[:idx] + 'cab' + s[idx+3:], True
    return s, False

def operation6(s):  # suffix 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
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

# Test the string
input_string = "bbcccbcbaccabaabaaa"
process_string(input_string)