def apply_rules(tokens):
    tokens = tokens.split()
    while True:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            pair = tokens[i] + " " + tokens[i+1]
            
            # Check if the current token ends with # and next token starts with #
            if tokens[i].endswith("#") and tokens[i+1].startswith("#"):
                if tokens[i][0] == "A" and tokens[i+1][-1] == "A":  # A# #A
                    tokens.pop(i)
                    tokens.pop(i)
                    changed = True
                    break
                elif tokens[i][0] == "B" and tokens[i+1][-1] == "B":  # B# #B
                    tokens.pop(i)
                    tokens.pop(i)
                    changed = True
                    break
                elif tokens[i][0] == "A" and tokens[i+1][-1] == "B":  # A# #B
                    tokens[i:i+2] = ["#B", "A#"]
                    changed = True
                    break
                elif tokens[i][0] == "B" and tokens[i+1][-1] == "A":  # B# #A
                    tokens[i:i+2] = ["#A", "B#"]
                    changed = True
                    break
            i += 1
            
        print(" ".join(tokens))
        if not changed:
            break
    
    return " ".join(tokens)

program = "#B B# #A B# A# A# #A #B B# B#"
print("Initial state:")
print(program)
print("\nTransformations:")
final_state = apply_rules(program)