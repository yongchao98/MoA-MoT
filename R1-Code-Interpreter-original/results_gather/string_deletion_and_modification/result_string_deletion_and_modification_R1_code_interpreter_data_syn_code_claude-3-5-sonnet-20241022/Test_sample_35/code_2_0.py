def check_and_apply(s):
    print(f"\nChecking string: {s}")
    
    # Operation 1: 'ca' not at start
    pos = s[1:].find('ca')
    print(f"Op1 - 'ca' after first char position: {pos}")
    if pos != -1:
        new_s = s[:pos+1] + s[pos+3:]
        print(f"Op1 applies: {s} -> {new_s}")
        return new_s, True
    
    # Operation 2: even number of b's
    b_count = s.count('b')
    print(f"Op2 - Count of 'b's: {b_count}")
    if b_count > 0 and b_count % 2 == 0:
        new_s = s + 'ab'
        print(f"Op2 applies: {s} -> {new_s}")
        return new_s, True
    
    # Operation 3: contains 'bca'
    pos = s.find('bca')
    print(f"Op3 - 'bca' position: {pos}")
    if pos != -1:
        new_s = s[:pos] + s[pos+3:]
        print(f"Op3 applies: {s} -> {new_s}")
        return new_s, True
    
    # Operation 4: prefix 'cb'
    print(f"Op4 - starts with 'cb': {s.startswith('cb')}")
    if s.startswith('cb'):
        new_s = 'aa' + s[2:-1]
        print(f"Op4 applies: {s} -> {new_s}")
        return new_s, True
    
    # Operation 5: prefix 'ab'
    print(f"Op5 - starts with 'ab': {s.startswith('ab')}")
    if s.startswith('ab'):
        new_s = 'ca' + s[2:]
        print(f"Op5 applies: {s} -> {new_s}")
        return new_s, True
    
    # Operation 6: suffix 'cc'
    print(f"Op6 - ends with 'cc': {s.endswith('cc')}")
    if s.endswith('cc'):
        new_s = 'a' + s[:-2] + 'b'
        print(f"Op6 applies: {s} -> {new_s}")
        return new_s, True
    
    return s, False

s = "aaccbcbaacccbbab"
changed = True
step = 1

while changed:
    print(f"\nStep {step}:")
    s, changed = check_and_apply(s)
    if changed:
        step += 1

print("\nFinal string:", s)