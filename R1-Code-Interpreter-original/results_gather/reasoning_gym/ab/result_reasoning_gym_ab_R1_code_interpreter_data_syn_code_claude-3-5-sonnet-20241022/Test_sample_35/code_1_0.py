def apply_rules(tokens):
    for i in range(len(tokens) - 1):
        # Only consider pairs where # symbols face each other
        if tokens[i].endswith("#") and tokens[i+1].startswith("#"):
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
while True:
    new_program, changed = apply_rules(program)
    if not changed:
        break
    program = new_program
    print(" ".join(program))  # Print each step

print("\nFinal state:")
print(" ".join(program))