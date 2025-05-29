def compute_program(tokens):
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i + 1]
            
            if pair == "A# #A" or pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                i = max(i - 1, 0)  # Move back one step to recheck the new pair
            elif pair == "A# #B":
                tokens[i], tokens[i + 1] = "#B", "A#"
                changed = True
                i += 1
            elif pair == "B# #A":
                tokens[i], tokens[i + 1] = "#A", "B#"
                changed = True
                i += 1
            else:
                i += 1
        
        if not changed:
            break
    
    return tokens

# Initial program
program = "#A B# #A B# #B B# B# A# #B B#".split()

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")