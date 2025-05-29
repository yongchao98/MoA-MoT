def apply_rules(program):
    tokens = program.split()
    i = 0
    
    while i < len(tokens) - 1:
        # Get current pair
        pair = f"{tokens[i]} {tokens[i+1]}"
        
        # Apply rules
        if pair == "A# #A":
            # Remove both tokens
            tokens.pop(i)
            tokens.pop(i)
            # Move back one position if possible to catch new pairs
            i = max(0, i-1)
        elif pair == "A# #B":
            tokens[i] = "#B"
            tokens[i+1] = "A#"
            # Stay at current position
            continue
        elif pair == "B# #A":
            tokens[i] = "#A"
            tokens[i+1] = "B#"
            # Stay at current position
            continue
        elif pair == "B# #B":
            # Remove both tokens
            tokens.pop(i)
            tokens.pop(i)
            # Move back one position if possible to catch new pairs
            i = max(0, i-1)
        else:
            # Move to next position if no rule applies
            i += 1
    
    return " ".join(tokens)

# Initial program
program = "B# B# #B A# #A B# B# #A #A #A"
result = apply_rules(program)
print(result)