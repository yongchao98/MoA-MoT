def apply_rules(tokens):
    # Returns (new_tokens, changed) where changed indicates if any rule was applied
    
    for i in range(len(tokens) - 1):
        pair = tokens[i] + " " + tokens[i+1]
        
        # Rule 1: A# #A becomes nothing
        if pair == "A# #A":
            return tokens[:i] + tokens[i+2:], True
            
        # Rule 2: A# #B becomes #B A#
        elif pair == "A# #B":
            return tokens[:i] + ["#B", "A#"] + tokens[i+2:], True
            
        # Rule 3: B# #A becomes #A B#
        elif pair == "B# #A":
            return tokens[:i] + ["#A", "B#"] + tokens[i+2:], True
            
        # Rule 4: B# #B becomes nothing
        elif pair == "B# #B":
            return tokens[:i] + tokens[i+2:], True
            
    return tokens, False

# Initial program
program = ["A#", "A#", "#B", "#B", "#B", "#A", "#B", "#B", "#A", "B#"]

# Keep applying rules until no more changes
steps = 0
while True:
    new_program, changed = apply_rules(program)
    if not changed:
        break
    program = new_program
    steps += 1
    
print(" ".join(program))