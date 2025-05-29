def apply_operations(s):
    operations_performed = True
    steps = []
    steps.append(s)
    
    while operations_performed:
        operations_performed = False
        
        # Operation 1: Replace first 'acb' with 'bca'
        if 'acb' in s:
            idx = s.index('acb')
            s = s[:idx] + 'bca' + s[idx+3:]
            operations_performed = True
            steps.append(f"1: {s}")
            continue
            
        # Operation 2: Delete first 'bca'
        if 'bca' in s:
            idx = s.index('bca')
            s = s[:idx] + s[idx+3:]
            operations_performed = True
            steps.append(f"2: {s}")
            continue
            
        # Operation 3: Remove middle character if length > 15
        if len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
            operations_performed = True
            steps.append(f"3: {s}")
            continue
            
        # Operation 4: Replace 'ca' prefix with 'bb' and append 'c'
        if s.startswith('ca'):
            s = 'bb' + s[2:] + 'c'
            operations_performed = True
            steps.append(f"4: {s}")
            continue
            
        # Operation 5: Delete 'bb' suffix
        if s.endswith('bb'):
            s = s[:-2]
            operations_performed = True
            steps.append(f"5: {s}")
            continue
            
        # Operation 6: Delete 'bc' prefix and append 'aa'
        if s.startswith('bc'):
            s = s[2:] + 'aa'
            operations_performed = True
            steps.append(f"6: {s}")
            continue
    
    print("\n".join(steps))
    print("\nFinal string:", s)

# Test with the given string
initial_string = "bcaccbaaabaacacbb"
apply_operations(initial_string)