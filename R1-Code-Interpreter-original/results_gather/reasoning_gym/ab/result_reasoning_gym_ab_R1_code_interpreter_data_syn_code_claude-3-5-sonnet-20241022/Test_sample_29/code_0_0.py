def apply_rules(program):
    # Convert program string to list of tokens
    tokens = program.split()
    
    def try_apply_rule(i):
        if i >= len(tokens) - 1:
            return False
            
        pair = f"{tokens[i]} {tokens[i+1]}"
        
        if pair == "A# #A" or pair == "B# #B":
            # Remove both tokens
            tokens.pop(i)
            tokens.pop(i)
            return True
        elif pair == "A# #B":
            # Swap and transform
            tokens[i:i+2] = ["#B", "A#"]
            return True
        elif pair == "B# #A":
            # Swap and transform
            tokens[i:i+2] = ["#A", "B#"]
            return True
        return False

    changed = True
    while changed:
        changed = False
        i = 0
        while i < len(tokens):
            if try_apply_rule(i):
                changed = True
            else:
                i += 1
        
    return " ".join(tokens)

# Initial program
program = "#B B# B# #A #A #A A# #B #A #A"
final_state = apply_rules(program)
print(final_state)