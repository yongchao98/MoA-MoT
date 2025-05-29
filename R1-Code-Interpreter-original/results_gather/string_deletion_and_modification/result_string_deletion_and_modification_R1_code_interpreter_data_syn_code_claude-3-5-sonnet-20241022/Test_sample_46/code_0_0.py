def apply_rules(s):
    def rule1(s):  # Replace 'ac' at start with 'zz'
        return 'zz' + s[2:] if s.startswith('ac') else s
    
    def rule2(s):  # Delete first 'bca'
        return s.replace('bca', '', 1) if 'bca' in s else s
    
    def rule3(s):  # Replace first 'acb' with 'bca'
        return s.replace('acb', 'bca', 1) if 'acb' in s else s
    
    def rule4(s):  # Replace ending 'ba' with 'ab'
        return s[:-2] + 'ab' if s.endswith('ba') else s
    
    def rule5(s):  # Remove middle char if length > 15
        if len(s) > 15:
            mid = len(s) // 2
            return s[:mid] + s[mid+1:]
        return s
    
    def rule6(s):  # Remove second char if starts with 'bb'
        return s[0] + s[2:] if s.startswith('bb') else s

    prev_s = ""
    current_s = s
    
    while prev_s != current_s:
        prev_s = current_s
        
        # Apply each rule in sequence
        current_s = rule1(current_s)
        print(f"After rule 1: {current_s}")
        
        current_s = rule2(current_s)
        print(f"After rule 2: {current_s}")
        
        current_s = rule3(current_s)
        print(f"After rule 3: {current_s}")
        
        current_s = rule4(current_s)
        print(f"After rule 4: {current_s}")
        
        current_s = rule5(current_s)
        print(f"After rule 5: {current_s}")
        
        current_s = rule6(current_s)
        print(f"After rule 6: {current_s}\n")
    
    return current_s

initial_string = "acabaacbacba"
result = apply_rules(initial_string)
print(f"Final result: {result}")