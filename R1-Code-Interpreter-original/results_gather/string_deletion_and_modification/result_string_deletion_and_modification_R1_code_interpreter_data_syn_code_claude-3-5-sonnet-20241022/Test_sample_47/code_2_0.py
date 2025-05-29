def find_all_patterns(s):
    patterns = []
    # Check prefix patterns
    if s.startswith('aa'): patterns.append(('aa prefix', 1))
    if s.startswith('cb'): patterns.append(('cb prefix', 2))
    if s.startswith('cc'): patterns.append(('cc prefix', 3))
    
    # Check suffix patterns
    if s.endswith('bb'): patterns.append(('bb suffix', 4))
    if s.endswith('aa'): patterns.append(('aa suffix', 6))
    
    # Check 'acb' pattern
    pos = s.find('acb')
    if pos != -1: patterns.append(('acb at ' + str(pos), 5))
    
    return patterns

def apply_rule(s, rule_num):
    if rule_num == 1:  # aa prefix
        return s[1:]
    elif rule_num == 2:  # cb prefix
        return 'aa' + s[2:-1]
    elif rule_num == 3:  # cc prefix
        return s[2:]
    elif rule_num == 4:  # bb suffix
        return s[:-2]
    elif rule_num == 5:  # acb -> bca
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:]
    elif rule_num == 6:  # aa suffix -> cc
        return s[:-2] + 'cc'
    return s

def process_string(s):
    steps = [(s, 'Initial')]
    while True:
        patterns = find_all_patterns(s)
        if not patterns:
            break
            
        print(f"\nCurrent string: {s}")
        print("Available patterns:", patterns)
        
        # Apply first available rule
        pattern_desc, rule_num = patterns[0]
        new_s = apply_rule(s, rule_num)
        steps.append((new_s, f"Rule {rule_num}: {pattern_desc}"))
        s = new_s
    
    print("\nAll steps:")
    for i, (string, desc) in enumerate(steps):
        print(f"Step {i}: {string} ({desc})")
    return s

# Start from where we left off
s = "bcacabbabbabac"
final = process_string(s)
print(f"\nFinal result: \"{final}\"")