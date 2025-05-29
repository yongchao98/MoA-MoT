def apply_rules(tokens):
    tokens = tokens.split()
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair
            pair = tokens[i] + " " + tokens[i+1]
            
            # Check for reversed pairs too
            rev_pair = tokens[i+1] + " " + tokens[i]
            
            # Apply rules (including reversed patterns)
            if pair == "A# #A" or pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif pair == "A# #B":
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
                break
            elif pair == "B# #A":
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
                break
            # Add reversed pattern checks
            elif rev_pair == "B# #A":  # This catches #A B# pattern
                tokens[i:i+2] = ["B#", "#A"]
                changed = True
                break
            elif rev_pair == "A# #B":  # This catches #B A# pattern
                tokens[i:i+2] = ["A#", "#B"]
                changed = True
                break
            i += 1
            
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