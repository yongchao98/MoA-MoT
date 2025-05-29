def transform_string(s):
    def try_rules(s):
        # Rule 1: Replace first 'acb' with 'bca'
        if 'acb' in s:
            idx = s.index('acb')
            return s[:idx] + 'bca' + s[idx+3:]
        
        # Rule 2: Remove first char if prefix is 'aa'
        if s.startswith('aa'):
            return s[1:]
        
        # Rule 3: Replace 'ca' prefix with 'bb' and append 'c'
        if s.startswith('ca'):
            return 'bb' + s[2:] + 'c'
        
        # Rule 4: Remove second char if starts with 'bb'
        if s.startswith('bb'):
            return s[0] + s[2:]
        
        # Rule 5: Delete last two chars if ends with 'bb'
        if s.endswith('bb'):
            return s[:-2]
        
        # Rule 6: Delete first two chars if starts with 'bc' and append 'aa'
        if s.startswith('bc'):
            return s[2:] + 'aa'
        
        return None

    steps = [s]
    while True:
        new_s = try_rules(s)
        if new_s is None:
            break
        s = new_s
        steps.append(s)
    
    print("Transformation steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "ccacbaabacbcacbccc"
result = transform_string(initial)
print(f"\nFinal result: {result}")