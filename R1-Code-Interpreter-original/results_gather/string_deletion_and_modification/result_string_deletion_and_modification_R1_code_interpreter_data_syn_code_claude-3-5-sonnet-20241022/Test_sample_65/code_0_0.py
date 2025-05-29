def apply_rules(s):
    def rule1(s):  # 'ca' -> 'bb' + 'c'
        if s.startswith('ca'):
            return 'bb' + s[2:] + 'c'
        return s

    def rule2(s):  # 'ab' -> 'ca'
        if s.startswith('ab'):
            return 'ca' + s[2:]
        return s

    def rule3(s):  # 'aa' at end -> 'cc'
        if s.endswith('aa'):
            return s[:-2] + 'cc'
        return s

    def rule4(s):  # 'acb' -> 'bca'
        if 'acb' in s:
            i = s.index('acb')
            return s[:i] + 'bca' + s[i+3:]
        return s

    def rule5(s):  # length > 15 -> remove middle char
        if len(s) > 15:
            mid = len(s) // 2
            return s[:mid] + s[mid+1:]
        return s

    def rule6(s):  # 'abc' -> 'cab'
        if 'abc' in s:
            i = s.index('abc')
            return s[:i] + 'cab' + s[i+3:]
        return s

    original = s
    s = rule1(s)
    if s != original: 
        print(f"Rule 1: {original} -> {s}")
        return s
    s = rule2(s)
    if s != original:
        print(f"Rule 2: {original} -> {s}")
        return s
    s = rule3(s)
    if s != original:
        print(f"Rule 3: {original} -> {s}")
        return s
    s = rule4(s)
    if s != original:
        print(f"Rule 4: {original} -> {s}")
        return s
    s = rule5(s)
    if s != original:
        print(f"Rule 5: {original} -> {s}")
        return s
    s = rule6(s)
    if s != original:
        print(f"Rule 6: {original} -> {s}")
        return s
    return s

# Initial string
s = "bacacbccabbcaabacaab"
prev = ""
while s != prev:
    prev = s
    s = apply_rules(s)

print("\nFinal string:", s)