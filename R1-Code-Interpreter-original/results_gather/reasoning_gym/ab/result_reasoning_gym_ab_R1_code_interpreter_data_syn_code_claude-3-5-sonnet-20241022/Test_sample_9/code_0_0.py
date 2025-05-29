def apply_rules(program):
    # Convert string to list for easier manipulation
    tokens = program.split()
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Check each pair of adjacent tokens
            current = tokens[i]
            next_token = tokens[i + 1]
            
            # Rule 1: A# #A -> nothing
            if current == "A#" and next_token == "#A":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            
            # Rule 2: A# #B -> #B A#
            elif current == "A#" and next_token == "#B":
                tokens[i] = "#B"
                tokens[i + 1] = "A#"
                changed = True
                break
            
            # Rule 3: B# #A -> #A B#
            elif current == "B#" and next_token == "#A":
                tokens[i] = "#A"
                tokens[i + 1] = "B#"
                changed = True
                break
            
            # Rule 4: B# #B -> nothing
            elif current == "B#" and next_token == "#B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            
            i += 1
        
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "#B #A B# B# #B #A #A #B A# B#"
final_state = apply_rules(program)
print(final_state)