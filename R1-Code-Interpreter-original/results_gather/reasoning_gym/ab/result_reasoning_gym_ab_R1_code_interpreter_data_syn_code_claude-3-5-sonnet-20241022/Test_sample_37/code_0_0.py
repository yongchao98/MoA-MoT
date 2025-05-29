def apply_rules(tokens):
    # Returns (new_tokens, changed) where changed indicates if any rule was applied
    i = 0
    while i < len(tokens) - 1:
        # Get current pair of tokens
        pair = tokens[i:i+2]
        
        # Check each rule
        if pair == ['A#', '#A']:
            # Remove both tokens
            tokens = tokens[:i] + tokens[i+2:]
            return tokens, True
        elif pair == ['A#', '#B']:
            # Swap and move '#' positions
            tokens = tokens[:i] + ['#B', 'A#'] + tokens[i+2:]
            return tokens, True
        elif pair == ['B#', '#A']:
            # Swap and move '#' positions
            tokens = tokens[:i] + ['#A', 'B#'] + tokens[i+2:]
            return tokens, True
        elif pair == ['B#', '#B']:
            # Remove both tokens
            tokens = tokens[:i] + tokens[i+2:]
            return tokens, True
        i += 1
    return tokens, False

# Initial program
program = ['A#', 'A#', 'B#', 'A#', '#A', '#A', '#B', '#B', '#B', '#B']

# Keep applying rules until no changes occur
while True:
    print("Current state:", ' '.join(program))
    new_program, changed = apply_rules(program)
    if not changed:
        break
    program = new_program

print("\nFinal state:", ' '.join(program))