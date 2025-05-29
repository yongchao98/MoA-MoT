def apply_rules(tokens):
    # Convert list to string for easier printing
    def tokens_to_str(t):
        return ' '.join(t)
    
    # Print initial state
    print("Initial:", tokens_to_str(tokens))
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            
            if pair == "A# #A" or pair == "B# #B":
                # Remove both tokens
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                print("After removal:", tokens_to_str(tokens))
                break
            
            elif pair == "A# #B":
                # Replace with "#B A#"
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
                print("After swap:", tokens_to_str(tokens))
                break
                
            elif pair == "B# #A":
                # Replace with "#A B#"
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
                print("After swap:", tokens_to_str(tokens))
                break
                
            i += 1
            
        if not changed:
            break
    
    return tokens_to_str(tokens)

# Initial program
program = "#A #A #B B# A# #A A# A# B# #B".split()
final = apply_rules(program)
print("\nFinal state:", final)