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
                new_tokens = rules[pair]
                tokens = tokens[:i] + new_tokens + tokens[i + 2:]
                changed = True
            else:
                i += 1
    return tokens

initial_program = ['#A', 'A#', '#A', 'B#', 'B#', '#A', 'A#', 'B#', 'A#', 'B#']
final_state = compute_program(initial_program)
print(final_state)