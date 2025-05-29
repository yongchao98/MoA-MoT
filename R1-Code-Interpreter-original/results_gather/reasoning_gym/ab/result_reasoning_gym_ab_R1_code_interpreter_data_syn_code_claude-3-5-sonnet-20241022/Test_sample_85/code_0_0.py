def apply_rules(tokens):
    # Returns (new_tokens, changed)
    # Changed is True if any rule was applied
    
    for i in range(len(tokens)-1):
        pair = tokens[i] + " " + tokens[i+1]
        
        # Rule 1: A# #A -> nothing
        if pair == "A# #A":
            return tokens[:i] + tokens[i+2:], True
            
        # Rule 2: A# #B -> #B A#
        if pair == "A# #B":
            return tokens[:i] + ["#B", "A#"] + tokens[i+2:], True
            
        # Rule 3: B# #A -> #A B#
        if pair == "B# #A":
            return tokens[:i] + ["#A", "B#"] + tokens[i+2:], True
            
        # Rule 4: B# #B -> nothing
        if pair == "B# #B":
            return tokens[:i] + tokens[i+2:], True
            
    return tokens, False

# Initial program
program = ["B#", "B#", "B#", "B#", "#B", "B#", "#B", "A#", "#A", "#A"]

# Keep applying rules until no changes
steps = 0
while True:
    print(f"Step {steps}: {' '.join(program)}")
    new_program, changed = apply_rules(program)
    if not changed:
        break
    program = new_program
    steps += 1

print("\nFinal state:", ' '.join(program))