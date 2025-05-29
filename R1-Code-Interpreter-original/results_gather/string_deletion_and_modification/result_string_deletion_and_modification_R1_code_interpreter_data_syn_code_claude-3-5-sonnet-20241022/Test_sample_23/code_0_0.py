def apply_rules(s):
    original = s
    steps = []
    steps.append(f"Start: {s}")
    
    while True:
        changed = False
        
        # Rule 1: Delete first 'bca'
        if 'bca' in s:
            pos = s.find('bca')
            s = s[:pos] + s[pos+3:]
            steps.append(f"Rule 1 (remove bca): {s}")
            changed = True
            continue
            
        # Rule 2: Remove 'cc' from start
        if s.startswith('cc'):
            s = s[2:]
            steps.append(f"Rule 2 (remove cc from start): {s}")
            changed = True
            continue
            
        # Rule 3: Remove 'ca' after first char
        if len(s) > 1:
            pos = s[1:].find('ca')
            if pos != -1:
                pos += 1  # Adjust position to account for full string
                s = s[:pos] + s[pos+2:]
                steps.append(f"Rule 3 (remove ca after first char): {s}")
                changed = True
                continue
                
        # Rule 4: Replace ending 'ba' with 'ab'
        if s.endswith('ba'):
            s = s[:-2] + 'ab'
            steps.append(f"Rule 4 (replace ending ba with ab): {s}")
            changed = True
            continue
            
        # Rule 5: Delete ending 'bb'
        if s.endswith('bb'):
            s = s[:-2]
            steps.append(f"Rule 5 (remove ending bb): {s}")
            changed = True
            continue
            
        # Rule 6: Replace starting 'ab' with 'ca'
        if s.startswith('ab'):
            s = 'ca' + s[2:]
            steps.append(f"Rule 6 (replace starting ab with ca): {s}")
            changed = True
            continue
            
        if not changed:
            break
            
    print("\n".join(steps))
    print(f"\nFinal result: {s}")

# Test with the given string
apply_rules("babbcacacabba")