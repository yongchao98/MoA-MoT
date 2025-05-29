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
                # Swap and transform
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
                break
            elif pair == "B# #A":
                # Swap and transform
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
                break
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

# Initial program
program = "#B #A A# A# #A B# B# #B A# B#"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")