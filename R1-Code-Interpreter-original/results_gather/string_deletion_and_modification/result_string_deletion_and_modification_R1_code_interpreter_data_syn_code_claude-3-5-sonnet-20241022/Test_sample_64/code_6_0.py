def apply_rules(s):
    step = 1
    while True:
        print(f"Step {step}: {s}")
        original = s
        
        # Rule 1: suffix 'cc' -> prepend 'a' and replace with 'b'
        if s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
            print(f"Applied Rule 1: {original} -> {s}")
            step += 1
            continue
            
        # Rule 2: length > 15 -> remove middle character
        if len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
            print(f"Applied Rule 2: {original} -> {s}")
            step += 1
            continue
            
        # Rule 3: 'acb' -> 'bca'
        if 'acb' in s:
            pos = s.find('acb')
            s = s[:pos] + 'bca' + s[pos+3:]
            print(f"Applied Rule 3: {original} -> {s}")
            step += 1
            continue
            
        # Rule 4: prefix 'cb' -> 'aa' and remove last char
        if s.startswith('cb'):
            s = 'aa' + s[2:-1]
            print(f"Applied Rule 4: {original} -> {s}")
            step += 1
            continue
            
        # Rule 5: starts with 'bb' -> remove second char
        if s.startswith('bb'):
            s = s[0] + s[2:]
            print(f"Applied Rule 5: {original} -> {s}")
            step += 1
            continue
            
        # Rule 6: starts with 'aa' -> remove first char
        if s.startswith('aa'):
            s = s[1:]
            print(f"Applied Rule 6: {original} -> {s}")
            step += 1
            continue
            
        if s == original:  # No rules applied
            break
            
    return s

s = "baccabcaabbcbca"
result = apply_rules(s)
print(f"\nFinal result: {result}")