def operation1(s):
    # If starts with 'bb', remove second character
    if s.startswith('bb'):
        return s[0] + s[2:], True
    return s, False

def operation2(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def operation3(s):
    # If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    return s, False

def operation4(s):
    # If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa', True
    return s, False

def operation5(s):
    # If even number of 'b's (and at least one), append 'ab'
    count_b = s.count('b')
    if count_b > 0 and count_b % 2 == 0:
        return s + 'ab', True
    return s, False

def operation6(s):
    # If ends with 'bb', delete last two characters
    if s.endswith('bb'):
        return s[:-2], True
    return s, False

def apply_operations(s):
    operations = [operation1, operation2, operation3, operation4, operation5, operation6]
    steps = []
    steps.append(s)
    
    while True:
        changed = False
        for op in operations:
            new_s, was_changed = op(s)
            if was_changed:
                s = new_s
                steps.append(s)
                changed = True
                break
        if not changed:
            break
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    print("\nFinal string:", s)

# Initial string
initial = "ababcaccbbbccbaabcbc"
apply_operations(initial)