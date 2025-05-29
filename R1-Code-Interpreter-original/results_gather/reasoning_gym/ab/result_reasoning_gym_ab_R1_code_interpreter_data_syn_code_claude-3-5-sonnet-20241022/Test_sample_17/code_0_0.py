def apply_rules(program):
    # Convert string to list for easier manipulation
    tokens = program.split()
    
    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Check each pair of adjacent tokens
            pair = tokens[i] + " " + tokens[i+1]
            
            # Rule 1: A# #A -> nothing
            if pair == "A# #A":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
                
            # Rule 2: B# #B -> nothing
            elif pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
                
            # Rule 3: A# #B -> #B A#
            elif pair == "A# #B":
                tokens[i] = "#B"
                tokens[i+1] = "A#"
                changed = True
                break
                
            # Rule 4: B# #A -> #A B#
            elif pair == "B# #A":
                tokens[i] = "#A"
                tokens[i+1] = "B#"
                changed = True
                break
                
            i += 1
    
    return " ".join(tokens)

# Initial program
program = "A# #A #B A# #A A# A# #A B# #B"
result = apply_rules(program)
print(result)