def apply_rules(tokens):
    # Returns (new_tokens, changed)
    for i in range(len(tokens) - 1):
        pair = tokens[i] + " " + tokens[i+1]
        
        # Rule: A# #A -> nothing
        if pair == "A# #A":
            return tokens[:i] + tokens[i+2:], True
            
        # Rule: A# #B -> #B A#
        elif pair == "A# #B":
            return tokens[:i] + ["#B", "A#"] + tokens[i+2:], True
            
        # Rule: B# #A -> #A B#
        elif pair == "B# #A":
            return tokens[:i] + ["#A", "B#"] + tokens[i+2:], True
            
        # Rule: B# #B -> nothing
        elif pair == "B# #B":
            return tokens[:i] + tokens[i+2:], True
            
    return tokens, False

# Initial program
program = ["B#", "B#", "B#", "B#", "#A", "#B", "#B", "#A", "#A", "#B"]

# Keep applying rules until no changes
while True:
    program, changed = apply_rules(program)
    if not changed:
        break

print(" ".join(program))