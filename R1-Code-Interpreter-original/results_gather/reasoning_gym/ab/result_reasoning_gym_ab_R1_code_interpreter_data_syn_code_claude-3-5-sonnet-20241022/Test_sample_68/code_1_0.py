def process_program(program):
    tokens = program.split()
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            # Get current pair of tokens
            pair = tokens[i] + " " + tokens[i+1]
            
            # Apply rules
            if pair == "A# #A" or pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif pair == "A# #B":
                tokens[i], tokens[i+1] = "#B", "A#"
                changed = True
                break
            elif pair == "B# #A" or pair == "#B #A":  # Added the reversed case
                tokens[i], tokens[i+1] = "#A", "B#"
                changed = True
                break
            i += 1
            
        if not changed:
            break
    
    return " ".join(tokens)

program = "B# #B #B B# #B #A #B A# A# A#"
print(process_program(program))