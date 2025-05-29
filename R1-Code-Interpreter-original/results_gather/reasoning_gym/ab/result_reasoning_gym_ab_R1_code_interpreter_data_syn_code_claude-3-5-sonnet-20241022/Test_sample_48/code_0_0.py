def apply_rules(tokens):
    # Returns (new_tokens, changed)
    # where changed indicates if any rule was applied
    i = 0
    while i < len(tokens) - 1:
        pair = tokens[i] + " " + tokens[i+1]
        
        if pair == "A# #A" or pair == "B# #B":
            # Remove both tokens
            tokens = tokens[:i] + tokens[i+2:]
            return tokens, True
            
        elif pair == "A# #B":
            # Replace with "#B A#"
            tokens = tokens[:i] + ["#B", "A#"] + tokens[i+2:]
            return tokens, True
            
        elif pair == "B# #A":
            # Replace with "#A B#"
            tokens = tokens[:i] + ["#A", "B#"] + tokens[i+2:]
            return tokens, True
            
        i += 1
    
    return tokens, False

# Initial sequence
sequence = ["#A", "#B", "#B", "#A", "#A", "#B", "#A", "A#", "#B", "B#"]

# Keep applying rules until no changes occur
while True:
    sequence, changed = apply_rules(sequence)
    if not changed:
        break

print(" ".join(sequence))