def apply_rules(program):
    tokens = program.split()
    
    def try_apply_rule(i):
        if i >= len(tokens) - 1:
            return False
            
        pair = f"{tokens[i]} {tokens[i+1]}"
        
        # Only apply rules when # faces #
        if tokens[i].endswith('#') and tokens[i+1].startswith('#'):
            if pair == "A# #A" or pair == "B# #B":
                tokens.pop(i)
                tokens.pop(i)
                return True
            elif pair == "A# #B":
                tokens[i:i+2] = ["#B", "A#"]
                return True
            elif pair == "B# #A":
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

program = "B# #B #B #A #B A# A# B# #A B#"
print(f"Initial: {program}")
print(f"Final: {apply_rules(program)}")