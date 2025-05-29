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
            pair = (tokens[i], tokens[i + 1])
            if pair in rules:
                # Apply the rule
                replacement = rules[pair]
                tokens[i:i + 2] = replacement
                changed = True
                # Move to the next pair after the replacement
                i += len(replacement) - 1
            else:
                i += 1
    
    # Join the tokens back into a string
    final_state = ' '.join(tokens)
    print(f"<<<{final_state}>>>")

# Initial program
program = "A# B# B# #A A# B# #A A# #A #B"
compute_program(program)