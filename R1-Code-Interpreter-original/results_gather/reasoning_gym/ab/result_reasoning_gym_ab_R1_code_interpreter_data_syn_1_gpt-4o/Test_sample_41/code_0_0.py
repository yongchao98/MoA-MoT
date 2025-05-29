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
                tokens[i:i + 2] = replacement
                changed = True
                if replacement:
                    i += len(replacement) - 1
                else:
                    i -= 1
            i += 1
    return tokens

initial_program = ['#A', 'A#', '#B', 'B#', '#B', '#A', 'A#', '#B', '#B', 'B#']
final_state = compute_program(initial_program)
print(final_state)