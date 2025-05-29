def rule1(s):
    if 'abc' in s:
        idx = s.index('abc')
        return s[:idx] + 'cab' + s[idx+3:], True
    return s, False

def rule2(s):
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', True
    return s, False

def rule3(s):
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def rule4(s):
    if s.endswith('aa'):
        return s[:-2] + 'cc', True
    return s, False

def rule5(s):
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab', True
    return s, False

def rule6(s):
    if len(s) > 1:
        idx = s[1:].find('ca')
        if idx != -1:
            idx += 1  # Adjust for the slice
            return s[:idx] + s[idx+2:], True
    return s, False

def process_string(s):
    print(f"Starting with: {s}")
    step = 1
    while True:
        original = s
        
        # Try each rule in order
        for rule_num, rule in enumerate([rule1, rule2, rule3, rule4, rule5, rule6], 1):
            s, changed = rule(s)
            if changed:
                print(f"Step {step}: Rule {rule_num} applied -> {s}")
                step += 1
                break
        
        # If no rules were applied, we're done
        if s == original:
            break
    
    return s

# Process the input string
input_str = "bbabacbcaaccaaaaccba"
final_str = process_string(input_str)
print(f"\nFinal string: {final_str}")