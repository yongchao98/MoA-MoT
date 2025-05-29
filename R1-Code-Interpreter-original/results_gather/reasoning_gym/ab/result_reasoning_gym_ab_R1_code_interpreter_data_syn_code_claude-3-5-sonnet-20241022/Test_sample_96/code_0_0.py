def apply_rules(tokens):
    # Convert to list for easier manipulation
    tokens = tokens.split()
    
    # Keep applying rules until no changes can be made
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = f"{tokens[i]} {tokens[i+1]}"
            
            if pair == "A# #A" or pair == "B# #B":
                # Remove both tokens
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
            
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "B# #B #A B# #B A# B# #A #A #A"
print(f"Initial: {program}")

# Process step by step
result = apply_rules(program)
print(f"Final: {result}")