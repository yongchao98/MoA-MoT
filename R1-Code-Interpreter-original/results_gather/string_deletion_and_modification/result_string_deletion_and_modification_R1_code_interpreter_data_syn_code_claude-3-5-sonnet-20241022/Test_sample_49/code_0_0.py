def apply_operations(s):
    def operation1(s):  # If starts with 'bb', remove second character
        if s.startswith('bb'):
            return s[0] + s[2:], True
        return s, False

    def operation2(s):  # If starts with 'cc', remove first two characters
        if s.startswith('cc'):
            return s[2:], True
        return s, False

    def operation3(s):  # If ends with 'bb', delete last two characters
        if s.endswith('bb'):
            return s[:-2], True
        return s, False

    def operation4(s):  # If contains 'ca' (not at start), remove first occurrence
        if len(s) > 1 and 'ca' in s[1:]:
            idx = s[1:].index('ca') + 1
            return s[:idx] + s[idx+2:], True
        return s, False

    def operation5(s):  # If starts with 'ca', replace with 'bb' and append 'c'
        if s.startswith('ca'):
            return 'bb' + s[2:] + 'c', True
        return s, False

    def operation6(s):  # If ends with 'aa', replace with 'cc'
        if s.endswith('aa'):
            return s[:-2] + 'cc', True
        return s, False

    current = s
    print(f"Starting string: {current}")
    
    while True:
        changed = False
        
        # Try each operation in order
        operations = [operation1, operation2, operation3, operation4, operation5, operation6]
        for i, op in enumerate(operations, 1):
            new_s, was_changed = op(current)
            if was_changed:
                print(f"Operation {i}: {current} -> {new_s}")
                current = new_s
                changed = True
                break
        
        if not changed:
            break
    
    return current

# Test with the given string
input_string = "cabbabbbabbbaabc"
result = apply_operations(input_string)
print(f"\nFinal result: {result}")