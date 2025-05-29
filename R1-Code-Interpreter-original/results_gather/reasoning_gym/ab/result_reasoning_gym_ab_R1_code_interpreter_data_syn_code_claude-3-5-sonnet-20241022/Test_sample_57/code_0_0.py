def apply_rules(program):
    # Convert program string to list of tokens
    tokens = program.split()
    
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
program = "A# #A B# #A A# B# #B #A B# #A"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")