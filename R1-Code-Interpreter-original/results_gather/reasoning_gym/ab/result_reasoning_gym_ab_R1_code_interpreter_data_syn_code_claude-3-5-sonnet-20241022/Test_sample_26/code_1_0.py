def apply_rule(tokens):
    # Print current state for debugging
    print("Current state:", tokens)
    
    for i in range(len(tokens) - 1):
        t1, t2 = tokens[i], tokens[i + 1]
        
        # Only process when '#' symbols face each other
        if t1.endswith('#') and t2.startswith('#'):
            new_tokens = tokens[:i]
            
            if (t1, t2) in [('A#', '#A'), ('B#', '#B')]:
                new_tokens.extend(tokens[i + 2:])
            elif t1 + t2 == 'A##B':
                new_tokens.extend(['#B', 'A#'])
                new_tokens.extend(tokens[i + 2:])
            elif t1 + t2 == 'B##A':
                new_tokens.extend(['#A', 'B#'])
                new_tokens.extend(tokens[i + 2:])
            
            return new_tokens, True
    return tokens, False

program = ['B#', 'B#', 'B#', 'A#', 'B#', '#B', '#B', '#A', '#A', '#A']
steps = 0

while True:
    program, changed = apply_rule(program)
    if not changed or steps > 20:  # Add step limit to prevent infinite loops
        break
    steps += 1

print("\nFinal state:", ' '.join(program))