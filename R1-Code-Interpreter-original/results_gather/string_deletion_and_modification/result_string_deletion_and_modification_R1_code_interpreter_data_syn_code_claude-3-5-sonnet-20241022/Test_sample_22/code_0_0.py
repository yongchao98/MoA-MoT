def apply_operations(s):
    def op1(s):  # length > 15, remove middle
        if len(s) > 15:
            mid = len(s) // 2
            return s[:mid] + s[mid+1:], True
        return s, False

    def op2(s):  # prefix 'cb', replace with 'aa' and delete last
        if s.startswith('cb'):
            return 'aa' + s[2:-1], True
        return s, False

    def op3(s):  # starts with 'bb', remove second char
        if s.startswith('bb'):
            return s[0] + s[2:], True
        return s, False

    def op4(s):  # prefix 'aa', remove first char
        if s.startswith('aa'):
            return s[1:], True
        return s, False

    def op5(s):  # contains 'ca' not at start, remove first occurrence
        if 'ca' in s[1:]:
            idx = s[1:].index('ca') + 1
            return s[:idx] + s[idx+2:], True
        return s, False

    def op6(s):  # prefix 'bc', delete first two and append 'aa'
        if s.startswith('bc'):
            return s[2:] + 'aa', True
        return s, False

    operations = [op1, op2, op3, op4, op5, op6]
    current = s
    steps = []
    steps.append(current)
    
    while True:
        changed = False
        for op in operations:
            new_s, was_changed = op(current)
            if was_changed:
                current = new_s
                steps.append(current)
                changed = True
                break
        if not changed:
            break
    
    print(f'Initial string: {s}')
    print('Steps:')
    for i, step in enumerate(steps):
        print(f'Step {i}: {step}')
    print(f'Final string: {current}')

# Test with the given string
input_string = "cabbbbbbbbbcaaaacb"
apply_operations(input_string)