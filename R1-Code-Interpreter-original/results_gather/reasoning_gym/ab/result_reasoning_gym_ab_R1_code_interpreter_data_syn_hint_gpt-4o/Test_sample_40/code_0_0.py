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
            i = max(i - 1, 0)  # Step back to recheck the previous pair
        else:
            i += 1
    
    return tokens

# Initial program
program = ['A#', 'B#', '#A', '#A', '#A', 'B#', '#A', 'B#', '#B', 'B#']

# Compute the final state
final_state = compute_program(program)

# Print the final state
print(final_state)