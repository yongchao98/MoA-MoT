def check_operations(s):
    print(f"Checking string: {s}")
    
    # Operation 1
    print(f"1. Ends with 'ca'? {s[-2:] if len(s) >= 2 else 'too short'}")
    
    # Operation 2
    print(f"2. Contains 'bca'? {s.find('bca') != -1}")
    
    # Operation 3
    print(f"3. Ends with 'ba'? {s[-2:] if len(s) >= 2 else 'too short'}")
    if s.endswith('ba'):
        new_s = s[:-2] + 'ab'
        print(f"   Can apply operation 3! Result would be: {new_s}")
    
    # Operation 4
    print(f"4. Starts with 'cb'? {s[:2] if len(s) >= 2 else 'too short'}")
    
    # Operation 5
    if len(s) > 1:
        ca_pos = s[1:].find('ca')
        print(f"5. Contains 'ca' after first char? {ca_pos != -1}")
    
    # Operation 6
    print(f"6. Starts with 'aa'? {s[:2] if len(s) >= 2 else 'too short'}")

check_operations("bacbcbacc")