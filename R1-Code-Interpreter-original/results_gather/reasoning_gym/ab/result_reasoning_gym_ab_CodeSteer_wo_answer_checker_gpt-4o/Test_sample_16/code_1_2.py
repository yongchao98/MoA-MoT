def compute_program(tokens):
    i = 0
    while i < len(tokens) - 1:
        pair = tokens[i] + " " + tokens[i + 1]
        print(f"Current tokens: {' '.join(tokens)}")
        print(f"Checking pair: {pair} at index {i}")
        
        if pair == "A# #A" or pair == "B# #B":
            print(f"Applying rule: {pair} -> nothing")
            tokens.pop(i)
            tokens.pop(i)
            i = max(i - 1, 0)  # Move back one step to recheck the new pair
        elif pair == "A# #B":
            print(f"Applying rule: {pair} -> #B A#")
            tokens[i], tokens[i + 1] = "#B", "A#"
            i += 1
        elif pair == "B# #A":
            print(f"Applying rule: {pair} -> #A B#")
            tokens[i], tokens[i + 1] = "#A", "B#"
            i += 1
        else:
            i += 1
        print(f"Tokens after operation: {' '.join(tokens)}\n")
    
    return tokens

# Initial program
program = "#A B# #A B# #B B# B# A# #B B#".split()

# Compute the final state
final_state = compute_program(program)

# Output the final state
print(f"<<<{' '.join(final_state)}>>>")