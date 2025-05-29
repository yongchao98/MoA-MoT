def apply_rules(s):
    def print_rule_application(rule_num, old_s, new_s):
        print(f"Rule {rule_num}: {old_s} -> {new_s}")
        
    changes_made = True
    current = s
    step = 1
    
    print(f"Initial: {current}")
    
    while changes_made:
        changes_made = False
        old = current
        
        # Rule 1: ends with 'ba' -> 'ab'
        if current.endswith('ba'):
            current = current[:-2] + 'ab'
            print_rule_application(1, old, current)
            continue
            
        # Rule 2: starts with 'aa' -> remove first 'a'
        if current.startswith('aa'):
            current = current[1:]
            print_rule_application(2, old, current)
            continue
            
        # Rule 3: ends with 'ac' -> 'cb'
        if current.endswith('ac'):
            current = current[:-2] + 'cb'
            print_rule_application(3, old, current)
            continue
            
        # Rule 4: starts with 'cc' -> remove 'cc'
        if current.startswith('cc'):
            current = current[2:]
            print_rule_application(4, old, current)
            continue
            
        # Rule 5: starts with 'bc' -> remove 'bc' and add 'aa'
        if current.startswith('bc'):
            current = current[2:] + 'aa'
            print_rule_application(5, old, current)
            continue
            
        # Rule 6: ends with 'cc' -> add 'a' at start and replace 'cc' with 'b'
        if current.endswith('cc'):
            current = 'a' + current[:-2] + 'b'
            print_rule_application(6, old, current)
            continue
            
        # If we get here and no changes were made, we're done
        if current == old:
            changes_made = False
        else:
            changes_made = True
            step += 1
            
    print(f"Final: {current}")
    return current

# Test with the string
s = "accbabcbabcbb"
result = apply_rules(s)