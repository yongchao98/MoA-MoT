def apply_rules(tokens):
    # Convert to list for easier manipulation
    tokens = tokens.split()
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair of tokens
            pair = tokens[i] + " " + tokens[i+1]
            
            # Apply rules
            if pair == "A# #A" or pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif pair == "A# #B":
                tokens[i] = "#B"
                tokens[i+1] = "A#"
                changed = True
                break
            elif pair == "B# #A":
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                changed = True
                break
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "#B #A #B #B A# #A B# #A B# #A"
print("Initial:", program)
print("Final:", apply_rules(program))