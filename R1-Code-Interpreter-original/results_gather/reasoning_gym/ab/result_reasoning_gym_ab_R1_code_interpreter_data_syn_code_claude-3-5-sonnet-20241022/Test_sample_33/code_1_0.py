def can_apply_rule(pair):
    # Check if '#' symbols face each other
    if len(pair) != 2:
        return False
    
    # For first token, '#' should be at end, for second token, '#' should be at start
    first_hash_at_end = pair[0].endswith('#')
    second_hash_at_start = pair[1].startswith('#')
    
    return not (first_hash_at_end or second_hash_at_start)

def apply_rules(tokens):
    # Returns (new_tokens, changed, position) where position indicates where rule was applied
    i = 0
    while i < len(tokens) - 1:
        pair = tokens[i:i+2]
        
        if can_apply_rule(pair):
            if pair == ['A#', '#A'] or pair == ['B#', '#B']:
                tokens = tokens[:i] + tokens[i+2:]
                return tokens, True, i
            elif pair == ['A#', '#B']:
                tokens = tokens[:i] + ['#B', 'A#'] + tokens[i+2:]
                return tokens, True, i
            elif pair == ['B#', '#A']:
                tokens = tokens[:i] + ['#A', 'B#'] + tokens[i+2:]
                return tokens, True, i
        i += 1
    return tokens, False, -1

# Initial program
program = ['#B', '#A', '#A', 'B#', '#B', '#B', '#B', '#A', '#B', 'B#']

# Keep applying rules and show each step
step = 0
print(f"Step {step}: {' '.join(program)}")

while True:
    program, changed, pos = apply_rules(program)
    if not changed:
        break
    step += 1
    print(f"Step {step}: {' '.join(program)}")