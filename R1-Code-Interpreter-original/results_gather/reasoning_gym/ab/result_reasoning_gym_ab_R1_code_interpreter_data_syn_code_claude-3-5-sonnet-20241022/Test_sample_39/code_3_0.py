def apply_rules(program):
    tokens = program.split()
    changed = True
    
    while changed:
        changed = False
        i = 0
        while i < len(tokens) - 1:
            current = tokens[i]
            next_token = tokens[i + 1]
            
            # Rule: A# #A -> nothing
            if current == "A#" and next_token == "#A":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            
            # Rule: B# #B -> nothing
            elif current == "B#" and next_token == "#B":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
            
            # Rule: A# #B -> #B A#
            elif current == "A#" and next_token == "#B":
                tokens[i] = "#B"
                tokens[i + 1] = "A#"
                changed = True
                break
            
            # Rule: B# #A -> #A B#
            elif current == "B#" and next_token == "#A":
                tokens[i] = "#A"
                tokens[i + 1] = "B#"
                changed = True
                break
                
            # Also check for B# B# pattern
            elif current == "B#" and next_token == "B#":
                tokens.pop(i)
                tokens.pop(i)
                changed = True
                break
                
            i += 1
            
    return " ".join(tokens)

program = "#A A# #A B# B# #A A# B# A# B#"
result = apply_rules(program)
print(result)