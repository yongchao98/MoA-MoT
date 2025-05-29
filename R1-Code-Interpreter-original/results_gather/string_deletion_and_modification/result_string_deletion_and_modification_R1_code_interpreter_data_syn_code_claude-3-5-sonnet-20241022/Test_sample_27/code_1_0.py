def apply_rules(s):
    print(f"Checking string: {s}")
    # Rule 1: prefix 'bc' -> delete first two and append 'aa'
    if s.startswith('bc'):
        result = s[2:] + 'aa'
        print(f"Rule 1 applied: {result}")
        return result, True
    
    # Rule 2: prefix 'aa' -> remove first character
    if s.startswith('aa'):
        result = s[1:]
        print(f"Rule 2 applied: {result}")
        return result, True
    
    # Rule 3: prefix 'ac' -> replace with 'zz'
    if s.startswith('ac'):
        result = 'zz' + s[2:]
        print(f"Rule 3 applied: {result}")
        return result, True
    
    # Rule 4: even number of 'b's -> append 'ab'
    if 'b' in s and s.count('b') % 2 == 0:
        result = s + 'ab'
        print(f"Rule 4 applied: {result}")
        return result, True
    
    # Rule 5: prefix 'ab' -> replace with 'ca'
    if s.startswith('ab'):
        result = 'ca' + s[2:]
        print(f"Rule 5 applied: {result}")
        return result, True
    
    # Rule 6: prefix 'cb' -> replace with 'aa' and delete last char
    if s.startswith('cb'):
        result = 'aa' + s[2:-1]
        print(f"Rule 6 applied: {result}")
        return result, True
    
    return s, False

s = "cacbbaabbcab"
while True:
    new_s, changed = apply_rules(s)
    if not changed:
        break
    s = new_s
print("\nFinal string:", s)