def apply_rules(tokens):
    # Returns (new_tokens, changed)
    i = 0
    while i < len(tokens) - 1:
        pair = tokens[i] + " " + tokens[i+1]
        
        # Print current state and pair being checked for debugging
        print(f"Current state: {' '.join(tokens)}")
        print(f"Checking pair: {pair} at position {i}")
        
        if pair == "A# #A" or pair == "B# #B":
            # Remove both tokens
            tokens = tokens[:i] + tokens[i+2:]
            print(f"Applied annihilation rule")
            return tokens, True
            
        elif pair == "A# #B":
            # Replace with "#B A#"
            tokens = tokens[:i] + ["#B", "A#"] + tokens[i+2:]
            print(f"Applied A# #B -> #B A# rule")
            return tokens, True
            
        elif pair == "B# #A":
            # Replace with "#A B#"
            tokens = tokens[:i] + ["#A", "B#"] + tokens[i+2:]
            print(f"Applied B# #A -> #A B# rule")
            return tokens, True
            
        i += 1
    
    return tokens, False

# Initial sequence
sequence = ["#A", "#B", "#B", "#A", "#A", "#B", "#A", "A#", "#B", "B#"]

# Keep applying rules until no changes occur
step = 1
while True:
    print(f"\nStep {step}:")
    sequence, changed = apply_rules(sequence)
    if not changed:
        break
    step += 1

print("\nFinal state:", " ".join(sequence))