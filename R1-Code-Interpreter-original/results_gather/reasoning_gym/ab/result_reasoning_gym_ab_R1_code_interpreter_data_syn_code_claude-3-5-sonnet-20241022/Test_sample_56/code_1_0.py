def apply_rules(tokens):
    tokens = tokens.split()
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            
            # Make sure to check for both orientations of B# #B
            if pair == "A# #A" or pair == "B# #B" or (tokens[i].endswith("#") and tokens[i+1].startswith("#") and tokens[i][0] == tokens[i+1][-1]):
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
            
        print(" ".join(tokens) if tokens else "(empty)")
        
        if not changed:
            break
    
    return " ".join(tokens) if tokens else "(empty)"

program = "#B B# #A B# A# A# #A #B B# B#"
print("Initial state:")
print(program)
print("\nTransformations:")
final_state = apply_rules(program)