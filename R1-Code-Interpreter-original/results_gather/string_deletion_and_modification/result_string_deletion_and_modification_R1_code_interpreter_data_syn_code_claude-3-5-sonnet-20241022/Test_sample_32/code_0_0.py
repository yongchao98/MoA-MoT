def apply_rules(s):
    # Returns (new_string, rule_applied)
    
    # Rule 1: 'acb' -> 'bca'
    if 'acb' in s:
        pos = s.index('acb')
        return s[:pos] + 'bca' + s[pos+3:], 1
    
    # Rule 2: 'ca' not at start -> remove
    if 'ca' in s[1:]:
        pos = s[1:].index('ca') + 1
        return s[:pos] + s[pos+2:], 2
    
    # Rule 3: prefix 'ca' -> 'bb' + append 'c'
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c', 3
    
    # Rule 4: 'bca' -> delete
    if 'bca' in s:
        pos = s.index('bca')
        return s[:pos] + s[pos+3:], 4
    
    # Rule 5: ends with 'ba' -> 'ab'
    if s.endswith('ba'):
        return s[:-2] + 'ab', 5
    
    # Rule 6: starts with 'bb' -> remove second char
    if s.startswith('bb'):
        return s[0] + s[2:], 6
    
    return s, 0

# Initial string
s = "cacccaacaaaaa"
iterations = 0
prev_s = ""

print("Starting string:", s)
while s != prev_s:
    prev_s = s
    s, rule = apply_rules(s)
    if rule > 0:
        print(f"After rule {rule}: {s}")
    iterations += 1

print("\nFinal string:", s)