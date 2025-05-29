def apply_rule1(s):
    # If starts with 'ac', replace with 'zz'
    if s.startswith('ac'):
        return 'zz' + s[2:], True
    return s, False

def apply_rule2(s):
    # If even number of b's, append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab', True
    return s, False

def apply_rule3(s):
    # If starts with 'cc', remove first two chars
    if s.startswith('cc'):
        return s[2:], True
    return s, False

def apply_rule4(s):
    # If length > 15, remove middle char
    if len(s) > 15:
        mid = len(s) // 2
        return s[:mid] + s[mid+1:], True
    return s, False

def apply_rule5(s):
    # If contains 'ca' (not at start), remove first occurrence after first char
    if 'ca' in s[1:]:
        idx = s.find('ca', 1)
        return s[:idx] + s[idx+2:], True
    return s, False

def apply_rule6(s):
    # If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        idx = s.find('acb')
        return s[:idx] + 'bca' + s[idx+3:], True
    return s, False

def process_string(s):
    print(f"Starting with: {s}")
    iteration = 1
    while True:
        original = s
        
        # Try each rule in order
        for rule_num, rule in enumerate([apply_rule1, apply_rule2, apply_rule3, 
                                       apply_rule4, apply_rule5, apply_rule6], 1):
            s, changed = rule(s)
            if changed:
                print(f"Iteration {iteration}, Rule {rule_num}: {s}")
                iteration += 1
                break
        
        # If no rules applied, we're done
        if s == original:
            break
    
    return s

# Process the input string
input_str = "babaabaaacaabaabb"
final_str = process_string(input_str)
print(f"\nFinal string: {final_str}")