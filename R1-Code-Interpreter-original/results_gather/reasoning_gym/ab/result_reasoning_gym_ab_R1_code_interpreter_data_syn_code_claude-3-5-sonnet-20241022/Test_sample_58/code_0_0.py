def apply_rules(tokens):
    def try_apply_rule(i):
        if i + 1 >= len(tokens):
            return False
        
        pair = tokens[i] + " " + tokens[i+1]
        
        if pair == "A# #A" or pair == "B# #B":
            # Remove both tokens
            tokens.pop(i)
            tokens.pop(i)
            return True
        elif pair == "A# #B":
            # Replace with "#B A#"
            tokens[i] = "#B"
            tokens[i+1] = "A#"
            return True
        elif pair == "B# #A":
            # Replace with "#A B#"
            tokens[i] = "#A"
            tokens[i+1] = "B#"
            return True
        return False

    # Initial sequence
    sequence = "B# A# #B #A #A #B #A #A #B #A"
    tokens = sequence.split()
    
    # Keep applying rules until no changes occur
    while True:
        change_made = False
        for i in range(len(tokens)-1):
            if try_apply_rule(i):
                change_made = True
                break
        if not change_made:
            break
    
    return " ".join(tokens)

print(apply_rules("B# A# #B #A #A #B #A #A #B #A"))