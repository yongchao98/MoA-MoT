def apply_rules(program):
    # Convert string to list for easier manipulation
    tokens = program.split()
    
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair of tokens
            pair = tokens[i] + " " + tokens[i+1]
            
            # Apply rules
            if pair == "A# #A" or pair == "B# #B":
                # Remove both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
            elif pair == "A# #B":
                # Swap tokens
                tokens[i], tokens[i+1] = "#B", "A#"
                changed = True
            elif pair == "B# #A":
                # Swap tokens
                tokens[i], tokens[i+1] = "#A", "B#"
                changed = True
            else:
                i += 1
    
    return " ".join(tokens)

# Initial program
program = "A# #B B# B# #B #B B# B# B# #B"
result = apply_rules(program)
print(result)