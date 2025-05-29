def rule1(s):
    # If 'ca' not at start, remove first occurrence after first char
    if len(s) > 2:
        pos = s[1:].find('ca')
        if pos != -1:
            return s[:pos+1] + s[pos+3:], True
    return s, False

def rule2(s):
    # Replace first 'abc' with 'cab'
    if 'abc' in s:
        return s.replace('abc', 'cab', 1), True
    return s, False

def rule3(s):
    # If suffix is 'bb', delete last two chars
    if s.endswith('bb'):
        return s[:-2], True
    return s, False

def rule4(s):
    # If suffix is 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    return s, False

def rule5(s):
    # If prefix is 'aa', remove first char
    if s.startswith('aa'):
        return s[1:], True
    return s, False

def rule6(s):
    # If starts with 'ac', replace first two chars with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:], True
    return s, False

def process_string(s):
    print(f"Starting with: {s}")
    rules = [rule1, rule2, rule3, rule4, rule5, rule6]
    step = 1
    
    while True:
        any_change = False
        for rule_num, rule in enumerate(rules, 1):
            new_s, changed = rule(s)
            if changed:
                print(f"Step {step}: Rule {rule_num} applied: {s} -> {new_s}")
                s = new_s
                any_change = True
                step += 1
                break
        if not any_change:
            break
    
    print(f"Final string: {s}")
    return s

# Process the given string
input_string = "aacbbacacabbba"
final_string = process_string(input_string)