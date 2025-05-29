def apply_rules(s):
    def print_step(step, rule, old, new):
        print(f"Step {step}: Rule {rule}: {old} -> {new}")
    
    step = 0
    print(f"Start: {s}")
    current = s
    
    while True:
        changed = False
        old = current
        
        # Check each rule in order, but look for patterns anywhere in the string
        
        # Rule 1: 'ba' -> 'ab'
        if 'ba' in current:
            pos = current.find('ba')
            current = current[:pos] + 'ab' + current[pos+2:]
            step += 1
            print_step(step, 1, old, current)
            continue
            
        # Rule 2: 'aa' prefix -> remove first 'a'
        if current.startswith('aa'):
            current = current[1:]
            step += 1
            print_step(step, 2, old, current)
            continue
            
        # Rule 3: 'ac' suffix -> 'cb'
        if current.endswith('ac'):
            current = current[:-2] + 'cb'
            step += 1
            print_step(step, 3, old, current)
            continue
            
        # Rule 4: 'cc' prefix -> remove 'cc'
        if current.startswith('cc'):
            current = current[2:]
            step += 1
            print_step(step, 4, old, current)
            continue
            
        # Rule 5: 'bc' anywhere -> remove and append 'aa'
        if 'bc' in current:
            pos = current.find('bc')
            current = current[:pos] + current[pos+2:] + 'aa'
            step += 1
            print_step(step, 5, old, current)
            continue
            
        # Rule 6: 'cc' suffix -> prepend 'a' and replace with 'b'
        if current.endswith('cc'):
            current = 'a' + current[:-2] + 'b'
            step += 1
            print_step(step, 6, old, current)
            continue
            
        if current == old:
            break
    
    print(f"Final: {current}")
    return current

# Starting with where we left off
s = "accbababbaaaa"
result = apply_rules(s)