def compute_program(tokens):
    rules = {
        ('A#', '#A'): [],
        ('A#', '#B'): ['#B', 'A#'],
        ('B#', '#A'): ['#A', 'B#'],
        ('B#', '#B'): []
    }
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = (tokens[i], tokens[i+1])
            if pair in rules:
                replacement = rules[pair]
                if replacement != tokens[i:i+2]:  # Check if replacement changes the sequence
                    tokens = tokens[:i] + replacement + tokens[i+2:]
                    changed = True
                    i += len(replacement) - 1
                else:
                    i += 1
            else:
                i += 1
        if not changed:
            break
    return tokens

initial_program = ['A#', '#A', '#B', 'B#', '#B', '#A', '#A', 'A#', '#B', 'B#']
final_state = compute_program(initial_program)
print(final_state)