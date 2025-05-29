def apply_rules(program):
    # Convert string to list for easier manipulation
    tokens = program.split()
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        for i in range(len(tokens) - 1):
            # Check each pair of adjacent tokens
            pair = f"{tokens[i]} {tokens[i+1]}"
            
            # Rule 1: A# #A -> nothing
            if pair == "A# #A":
                tokens = tokens[:i] + tokens[i+2:]
                changed = True
                break
                
            # Rule 2: A# #B -> #B A#
            elif pair == "A# #B":
                tokens = tokens[:i] + ["#B", "A#"] + tokens[i+2:]
                changed = True
                break
                
            # Rule 3: B# #A -> #A B#
            elif pair == "B# #A":
                tokens = tokens[:i] + ["#A", "B#"] + tokens[i+2:]
                changed = True
                break
                
            # Rule 4: B# #B -> nothing
            elif pair == "B# #B":
                tokens = tokens[:i] + tokens[i+2:]
                changed = True
                break
        
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "B# #B #A B# B# #A B# A# #B #B"
final_state = apply_rules(program)
print(final_state)