def apply_rules(s):
    print(f"\nStarting with: {s}")
    step = 1
    
    while True:
        print(f"\nStep {step}:")
        print(f"Current string: {s}")
        original = s
        
        # Rule 3: Replace 'acb' with 'bca'
        if 'acb' in s:
            pos = s.find('acb')
            print(f"Found 'acb' at position {pos}: {s[:pos]}|{s[pos:pos+3]}|{s[pos+3:]}")
            s = s[:pos] + 'bca' + s[pos+3:]
            print(f"After replacing 'acb' with 'bca': {s}")
            step += 1
            continue
        
        # If no changes were made, we're done
        if s == original:
            break
    
    return s

s = "baccabcaabbcbca"
result = apply_rules(s)
print(f"\nFinal result: {result}")