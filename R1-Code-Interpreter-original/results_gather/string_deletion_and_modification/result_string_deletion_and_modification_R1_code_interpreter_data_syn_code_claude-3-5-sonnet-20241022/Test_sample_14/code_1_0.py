def apply_rules(s):
    def highlight_all_patterns(s):
        print("\nAnalyzing patterns in:", s)
        # Check for consecutive 'abc'
        for i in range(len(s)-2):
            substr = s[i:i+3]
            if substr == 'abc':
                print(f"Found 'abc' at position {i}: {s[:i]}[{substr}]{s[i+3:]}")
    
    print(f"Starting with: {s}")
    while True:
        original = s
        highlight_all_patterns(s)
        
        # Apply rules in order
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
            print("Rule 1 applied: ends with 'ba' -> 'ab'")
        elif 'bca' in s:
            idx = s.find('bca')
            s = s[:idx] + s[idx+3:]
            print("Rule 2 applied: removed 'bca'")
        elif s.startswith('bb'):
            s = s[0] + s[2:]
            print("Rule 3 applied: removed second char from 'bb'")
        elif s.startswith('aa'):
            s = s[1:]
            print("Rule 4 applied: removed first 'a' from 'aa'")
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
            print("Rule 5 applied: 'cb' -> 'aa' and removed last char")
        elif 'abc' in s:
            idx = s.find('abc')
            s = s[:idx] + 'cab' + s[idx+3:]
            print(f"Rule 6 applied: replaced 'abc' with 'cab' at position {idx}")
        
        if s == original:
            break
            
        print(f"String is now: {s}")
    
    return s

initial = "aaaaacabacbcbc"
final = apply_rules(initial)
print(f"\nFinal result: {final}")