def apply_rules(s):
    def print_step(step, rule, old, new):
        print(f"Step {step}: Rule {rule}: {old} -> {new}")
    
    step = 0
    print(f"Start: {s}")
    current = s
    
    while True:
        old = current
        changed = False
        
        # Rule 1: If the string ends with 'ba', replace it with 'ab'
        if current.endswith('ba'):
            current = current[:-2] + 'ab'
            step += 1
            print_step(step, 1, old, current)
            continue
            
        # Rule 2: If the string prefix is 'aa', remove the first character
        if current.startswith('aa'):
            current = current[1:]
            step += 1
            print_step(step, 2, old, current)
            continue
            
        # Rule 3: If the string suffix is 'ac', replace it with 'cb'
        if current.endswith('ac'):
            current = current[:-2] + 'cb'
            step += 1
            print_step(step, 3, old, current)
            continue
            
        # Rule 4: If the string starts with 'cc', remove the first two characters
        if current.startswith('cc'):
            current = current[2:]
            step += 1
            print_step(step, 4, old, current)
            continue
            
        # Rule 5: If the string starts with 'bc', delete first two chars and append 'aa'
        if current.startswith('bc'):
            current = current[2:] + 'aa'
            step += 1
            print_step(step, 5, old, current)
            continue
            
        # Rule 6: If the string ends with 'cc', replace with 'b' and prepend 'a'
        if current.endswith('cc'):
            current = 'a' + current[:-2] + 'b'
            step += 1
            print_step(step, 6, old, current)
            continue
            
        if current == old:
            break
    
    print(f"Final: {current}")
    return current

# Starting with our string
s = "accbababbaaaa"
result = apply_rules(s)