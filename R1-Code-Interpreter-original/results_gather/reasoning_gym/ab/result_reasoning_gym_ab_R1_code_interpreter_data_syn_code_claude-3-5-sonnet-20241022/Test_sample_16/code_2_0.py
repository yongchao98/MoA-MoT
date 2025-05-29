def apply_rules(program):
    tokens = program.split()
    
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            
            # Rule 1: A# #A becomes nothing
            if pair == "A# #A":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            
            # Rule 2: A# #B becomes #B A#
            elif pair == "A# #B":
                tokens[i:i+2] = ["#B", "A#"]
                changed = True
                break
            
            # Rule 3: B# #A becomes #A B#
            elif pair == "B# #A":
                tokens[i:i+2] = ["#A", "B#"]
                changed = True
                break
            
            # Rule 4: B# #B becomes nothing
            elif pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            
            i += 1
            
        if not changed:
            break
        
        print(f"Step: {' '.join(tokens)}")  # Debug print
    
    return " ".join(tokens)

program = "#A #A B# B# A# B#"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")