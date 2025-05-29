def apply_rules(tokens):
    # Convert to list for easier manipulation
    tokens = tokens.split()
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair
            pair = tokens[i] + " " + tokens[i+1]
            
            # Apply rules
            if pair == "A# #A" or pair == "B# #B":
                # Delete both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif pair == "A# #B":
                # Replace with "#B A#"
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
                break
            elif pair == "B# #A":
                # Replace with "#A B#"
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
                break
            i += 1
            
        # Print current state
        print("Current state:", " ".join(tokens))
        
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "#B #A #A B# #B #A A# #A B# A#"
print("Initial state:", program)
print("\nProcessing...")
final_state = apply_rules(program)
print("\nFinal state:", final_state)