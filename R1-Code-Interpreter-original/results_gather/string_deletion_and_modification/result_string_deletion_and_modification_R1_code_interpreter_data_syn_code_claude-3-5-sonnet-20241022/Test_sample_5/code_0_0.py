def apply_rule1(s):
    # If ends with 'ba', replace with 'ab'
    if s.endswith('ba'):
        return s[:-2] + 'ab'
    return s

def apply_rule2(s):
    # If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc'
    return s

def apply_rule3(s):
    # If prefix is 'ca', replace with 'bb' and append 'c'
    if s.startswith('ca'):
        return 'bb' + s[2:] + 'c'
    return s

def apply_rule4(s):
    # If prefix is 'bc', delete first two chars and append 'aa'
    if s.startswith('bc'):
        return s[2:] + 'aa'
    return s

def apply_rule5(s):
    # If length > 15, remove middle character
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:]
    return s

def apply_rule6(s):
    # If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:]
    return s

def apply_rules(s):
    print(f"Starting with: {s}")
    prev = ""
    while prev != s:
        prev = s
        
        # Apply each rule in sequence
        for rule_num, rule in enumerate([apply_rule1, apply_rule2, apply_rule3, 
                                       apply_rule4, apply_rule5, apply_rule6], 1):
            new_s = rule(s)
            if new_s != s:
                print(f"Rule {rule_num}: {s} -> {new_s}")
                s = new_s
                break
    
    print(f"Final string: {s}")
    return s

# Test with the given string
initial = "abbbabbbacbabbbccaca"
result = apply_rules(initial)