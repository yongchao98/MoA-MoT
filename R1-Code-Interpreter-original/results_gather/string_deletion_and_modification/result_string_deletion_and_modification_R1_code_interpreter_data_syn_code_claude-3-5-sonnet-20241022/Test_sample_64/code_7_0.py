def print_with_highlight(s, pos, length, rule_num):
    before = s[:pos]
    highlight = s[pos:pos+length]
    after = s[pos+length:]
    print(f"Rule {rule_num} applies: {before}|{highlight}|{after}")

def apply_rules(s):
    step = 1
    while True:
        print(f"\nStep {step}: {s}")
        original = s
        changed = False
        
        # Rule 1: suffix 'cc'
        if s.endswith('cc'):
            print_with_highlight(s, len(s)-2, 2, 1)
            s = 'a' + s[:-2] + 'b'
            print(f"-> {s}")
            changed = True
        
        # Rule 2: length > 15
        elif len(s) > 15:
            mid = len(s) // 2
            print_with_highlight(s, mid, 1, 2)
            s = s[:mid] + s[mid+1:]
            print(f"-> {s}")
            changed = True
        
        # Rule 3: contains 'acb'
        elif 'acb' in s:
            pos = s.find('acb')
            print_with_highlight(s, pos, 3, 3)
            s = s[:pos] + 'bca' + s[pos+3:]
            print(f"-> {s}")
            changed = True
        
        # Rule 4: starts with 'cb'
        elif s.startswith('cb'):
            print_with_highlight(s, 0, 2, 4)
            s = 'aa' + s[2:-1]
            print(f"-> {s}")
            changed = True
        
        # Rule 5: starts with 'bb'
        elif s.startswith('bb'):
            print_with_highlight(s, 0, 2, 5)
            s = s[0] + s[2:]
            print(f"-> {s}")
            changed = True
        
        # Rule 6: starts with 'aa'
        elif s.startswith('aa'):
            print_with_highlight(s, 0, 2, 6)
            s = s[1:]
            print(f"-> {s}")
            changed = True
        
        if not changed:
            break
        step += 1
    
    return s

initial = "baccabcaabbcbca"
result = apply_rules(initial)
print(f"\nFinal result: {result}")