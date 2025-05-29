def transform_string(s):
    def check_and_apply(s):
        # Print current state and all possible conditions
        print(f"\nChecking string: {s}")
        print(f"1. Ends with 'ac'? {s.endswith('ac')}")
        print(f"2. Starts with 'cc'? {s.startswith('cc')}")
        print(f"3. Contains 'abc'? {'abc' in s}")
        print(f"4. Starts with 'bc'? {s.startswith('bc')}")
        print(f"5. Contains 'ca' after first char? {'ca' in s[1:]}")
        print(f"6. Ends with 'bb'? {s.endswith('bb')}")
        
        # Rule 1: If suffix is 'ac', replace with 'cb'
        if s.endswith('ac'):
            return (s[:-2] + 'cb', 1)
            
        # Rule 2: If starts with 'cc', remove first two characters
        if s.startswith('cc'):
            return (s[2:], 2)
            
        # Rule 3: If contains 'abc', replace first occurrence with 'cab'
        if 'abc' in s:
            return (s.replace('abc', 'cab', 1), 3)
            
        # Rule 4: If prefix is 'bc', delete first two and append 'aa'
        if s.startswith('bc'):
            return (s[2:] + 'aa', 4)
            
        # Rule 5: If contains 'ca' not at start, remove first occurrence after first char
        if 'ca' in s[1:]:
            idx = s.find('ca', 1)
            return (s[:idx] + s[idx+2:], 5)
            
        # Rule 6: If suffix is 'bb', delete last two characters
        if s.endswith('bb'):
            return (s[:-2], 6)
            
        return (s, 0)  # No rule applied

    current = s
    step = 0
    print(f"Step {step}: {current}")
    
    while True:
        new_s, rule = check_and_apply(current)
        if rule == 0:  # No rule could be applied
            break
        step += 1
        print(f"Step {step}: Applied Rule {rule}")
        print(f"  {current} -> {new_s}")
        current = new_s
    
    return current

# Initial string
s = "bbccbbbbacbabbba"
result = transform_string(s)
print(f"\nFinal result: {result}")