def compute_program(tokens):
    rules = {
        ('A#', '#A'): [],
        ('A#', '#B'): ['#B', 'A#'],
        ('B#', '#A'): ['#A', 'B#'],
        ('B#', '#B'): []
    }
    
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = (tokens[i], tokens[i + 1])
            if pair in rules:
                replacement = rules[pair]
                tokens = tokens[:i] + replacement + tokens[i + 2:]
                changed = True
                i = max(0, i - 1)  # Step back to check for new pairs
            else:
                i += 1
    return tokens

# Initial program
program = ['A#', 'B#', 'B#', '#B', '#A', 'A#', 'B#', 'A#', 'B#', 'A#']
final_state = compute_program(program)
print(final_state)