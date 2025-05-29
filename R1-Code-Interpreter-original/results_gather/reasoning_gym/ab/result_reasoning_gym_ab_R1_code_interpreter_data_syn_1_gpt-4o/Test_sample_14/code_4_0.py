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
                i = 0  # Reset index to start from the beginning after a change
            else:
                i += 1
    return tokens

initial_program = ['A#', '#B', 'A#', 'B#', 'B#', 'A#', 'A#', 'A#', '#A', 'A#']
final_state = compute_program(initial_program)
print(final_state)