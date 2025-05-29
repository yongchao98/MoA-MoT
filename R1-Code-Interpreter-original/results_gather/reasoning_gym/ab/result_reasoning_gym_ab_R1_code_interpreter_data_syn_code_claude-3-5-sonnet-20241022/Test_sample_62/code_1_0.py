def process_program(program):
    tokens = program.split()
    step = 0
    
    while True:
        changed = False
        i = 0
        print(f"Step {step}: {' '.join(tokens)}")
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            
            if pair == "A# #A" or pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            elif pair == "A# #B":
                tokens[i], tokens[i+1] = "#B", "A#"
                changed = True
                break
            elif pair == "B# #A":
                tokens[i], tokens[i+1] = "#A", "B#"
                changed = True
                break
            i += 1
        
        if not changed:
            break
        step += 1
    
    return " ".join(tokens)

program = "B# A# #A #A #B A# #B #A #A #B"
print("\nFinal result:", process_program(program))