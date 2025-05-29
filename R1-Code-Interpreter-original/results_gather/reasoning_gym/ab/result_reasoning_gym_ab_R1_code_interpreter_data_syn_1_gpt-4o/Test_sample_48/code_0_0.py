def compute_program(tokens):
    # Define the rewrite rules
    rules = {
        ('A#', '#A'): [],
        ('A#', '#B'): ['#B', 'A#'],
        ('B#', '#A'): ['#A', 'B#'],
        ('B#', '#B'): []
    }
    
    # Convert the input string to a list of tokens
    tokens = tokens.split()
    
    # Flag to check if any rewrite was applied
    rewrite_applied = True
    
    # Continue applying rules until no more rewrites can be applied
    while rewrite_applied:
        rewrite_applied = False
        i = 0
        while i < len(tokens) - 1:
            pair = (tokens[i], tokens[i + 1])
            if pair in rules:
                # Apply the rule
                new_tokens = rules[pair]
                tokens = tokens[:i] + new_tokens + tokens[i + 2:]
                rewrite_applied = True
                # Move to the next pair after the rewrite
                i += len(new_tokens) - 1
            else:
                i += 1
    
    # Return the final state of the program
    return ' '.join(tokens)

# Initial program
initial_program = "#A #B #B #A #A #B #A A# #B B#"

# Compute the final state
final_state = compute_program(initial_program)
print(final_state)