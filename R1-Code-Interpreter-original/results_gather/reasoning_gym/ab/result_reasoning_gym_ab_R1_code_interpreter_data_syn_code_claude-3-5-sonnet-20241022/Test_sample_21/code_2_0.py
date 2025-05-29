def apply_rules(program):
    tokens = program.split()
    print("\nStep by step:")
    print(tokens)
    
    def try_apply_rule(i):
        if i >= len(tokens) - 1:
            return False
            
        t1, t2 = tokens[i], tokens[i+1]
        
        # Only apply rules when # faces # (one ends with #, next starts with #)
        if t1.endswith('#') and t2.startswith('#'):
            if (t1 == "A#" and t2 == "#A") or (t1 == "B#" and t2 == "#B"):
                tokens.pop(i)
                tokens.pop(i)
                print(f"Remove pair {t1} {t2}: {tokens}")
                return True
            elif t1 == "A#" and t2 == "#B":
                tokens[i:i+2] = ["#B", "A#"]
                print(f"Swap {t1} {t2}: {tokens}")
                return True
            elif t1 == "B#" and t2 == "#A":
                tokens[i:i+2] = ["#A", "B#"]
                print(f"Swap {t1} {t2}: {tokens}")
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