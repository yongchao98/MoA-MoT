def apply_rules(program):
    tokens = program.split()
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair of tokens
            current = tokens[i]
            next_token = tokens[i+1]
            
            # Check if '#' symbols are facing each other
            if current.endswith('#') and next_token.startswith('#'):
                pair = current + " " + next_token
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
            i += 1
        
        if not changed:
            break
    
    return " ".join(tokens)

program = "#A #A B# B# A# B#"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")