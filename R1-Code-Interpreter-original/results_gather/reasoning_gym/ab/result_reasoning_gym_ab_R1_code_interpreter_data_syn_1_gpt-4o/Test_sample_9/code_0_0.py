def compute_program(tokens):
    rules = {
        ('A#', '#A'): [],
        ('A#', '#B'): ['#B', 'A#'],
        ('B#', '#A'): ['#A', 'B#'],
        ('B#', '#B'): []
    }
    
    i = 0
    while i < len(tokens) - 1:
        pair = (tokens[i], tokens[i + 1])
        if pair in rules:
            replacement = rules[pair]
            tokens = tokens[:i] + replacement + tokens[i + 2:]
            i = max(i - 1, 0)  # Step back to check for new pairs formed
        else:
            i += 1
    
    return tokens

# Initial program
program = ['#B', '#A', 'B#', 'B#', '#B', '#A', '#A', '#B', 'A#', 'B#']
final_state = compute_program(program)
print(final_state)