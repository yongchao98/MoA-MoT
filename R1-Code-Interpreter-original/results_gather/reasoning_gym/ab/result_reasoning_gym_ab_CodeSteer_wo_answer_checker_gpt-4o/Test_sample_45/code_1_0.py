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
    
    # Keep applying rules until no changes occur
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = (tokens[i], tokens[i+1])
            if pair in rules:
                # Apply the rule
                new_tokens = rules[pair]
                tokens = tokens[:i] + new_tokens + tokens[i+2:]
                changed = True
                # Move to the next pair after the current one
                i += len(new_tokens) - 1
            else:
                i += 1
    
    # Join the tokens back into a string
    final_state = ' '.join(tokens)
    return final_state

# Initial program
program = "#A #A B# B# B# A# #A #B B# B#"

# Compute the final state
final_state = compute_program(program)

# Output the final state in the specified format
print(f"<<<{final_state}>>>")