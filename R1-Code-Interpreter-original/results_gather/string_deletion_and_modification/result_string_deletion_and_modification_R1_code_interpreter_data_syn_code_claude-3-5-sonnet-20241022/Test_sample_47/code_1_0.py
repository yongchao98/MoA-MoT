def apply_rules(s):
    # Returns (new_string, changed)
    
    # Rule 1: If prefix is 'aa', remove first character
    if s.startswith('aa'):
        return s[1:], True, "Rule 1"
    
    # Rule 2: If prefix is 'cb', replace with 'aa' and delete last character
    if s.startswith('cb'):
        return 'aa' + s[2:-1], True, "Rule 2"
    
    # Rule 3: If starts with 'cc', remove first two characters
    if s.startswith('cc'):
        return s[2:], True, "Rule 3"
    
    # Rule 4: If suffix is 'bb', delete last two characters
    if s.endswith('bb'):
        return s[:-2], True, "Rule 4"
    
    # Rule 5: If contains 'acb', replace first occurrence with 'bca'
    if 'acb' in s:
        pos = s.find('acb')
        return s[:pos] + 'bca' + s[pos+3:], True, "Rule 5"
    
    # Rule 6: If ends with 'aa', replace with 'cc'
    if s.endswith('aa'):
        return s[:-2] + 'cc', True, "Rule 6"
    
    return s, False, "No rule"

def process_string(s):
    steps = [s]
    rules_applied = []
    iteration = 0
    max_iterations = 100  # Safety limit
    
    while iteration < max_iterations:
        new_s, changed, rule = apply_rules(s)
        if not changed:
            break
        s = new_s
        steps.append(s)
        rules_applied.append(rule)
        iteration += 1
    
    print("Steps:")
    for i, (step, rule) in enumerate(zip(steps, ['Initial'] + rules_applied)):
        print(f"Step {i}: {step} ({rule})")
    return s

# Process the given string
initial = "bcacabbabbabac"  # Starting from where we left off
final = process_string(initial)
print(f"\nFinal result: \"{final}\"")