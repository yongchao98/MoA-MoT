def show_matches(s):
    print(f"\nAnalyzing string: {s}")
    print(f"Length: {len(s)}")
    
    # Check for 'acb'
    pos = 0
    while pos < len(s)-2:
        if s[pos:pos+3] == 'acb':
            print(f"Found 'acb' at position {pos}: {s[:pos]}|{s[pos:pos+3]}|{s[pos+3:]}")
        pos += 1
    
    # Check for other patterns
    if s.endswith('cc'):
        print(f"Ends with 'cc': {s[:-2]}|{s[-2:]}")
    if s.startswith('cb'):
        print(f"Starts with 'cb': {s[:2]}|{s[2:]}")
    if s.startswith('bb'):
        print(f"Starts with 'bb': {s[:2]}|{s[2:]}")
    if s.startswith('aa'):
        print(f"Starts with 'aa': {s[:2]}|{s[2:]}")

s = "baccabcaabbcbca"
show_matches(s)