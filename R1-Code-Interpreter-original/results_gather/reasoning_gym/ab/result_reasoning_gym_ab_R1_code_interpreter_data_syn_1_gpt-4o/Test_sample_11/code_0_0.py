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
    
    # Flag to check if any rule was applied
    changed = True
    
    # Continue applying rules until no changes occur
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = (tokens[i], tokens[i + 1])
            if pair in rules:
                # Apply the rule
                new_tokens = rules[pair]
                tokens = tokens[:i] + new_tokens + tokens[i + 2:]
                changed = True
                # Move to the next pair after the current one
                i += len(new_tokens) - 1
            else:
                i += 1
    
    return tokens

# Initial program
program = "A# #A B# #B #A A# #B #B A# #B"

# Compute the final state of the program
final_state = compute_program(program)

# Print the final state
print(" ".join(final_state))