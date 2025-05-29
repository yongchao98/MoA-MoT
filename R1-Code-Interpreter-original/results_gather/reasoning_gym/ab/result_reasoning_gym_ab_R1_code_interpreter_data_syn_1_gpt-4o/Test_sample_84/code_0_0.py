def compute_program(tokens):
    rules = {
        ('A#', '#A'): [],
        ('A#', '#B'): ['#B', 'A#'],
        ('B#', '#A'): ['#A', 'B#'],
        ('B#', '#B'): []
    }
    
    i = 0
    while i < len(tokens) - 1:
        pair = (tokens[i], tokens[i+1])
        if pair in rules:
            # Apply the rule
            replacement = rules[pair]
            tokens = tokens[:i] + replacement + tokens[i+2:]
            # Restart from the beginning
            i = 0
        else:
            i += 1
    
    return tokens

# Initial program
program = ['#A', 'A#', 'A#', 'B#', '#B', '#B', 'B#', '#A', '#B', '#B']

# Compute the final state
final_state = compute_program(program)

# Print the final state
print(final_state)