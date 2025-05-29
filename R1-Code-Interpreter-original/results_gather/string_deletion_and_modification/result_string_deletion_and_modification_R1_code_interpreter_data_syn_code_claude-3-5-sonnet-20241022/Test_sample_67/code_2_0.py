def find_acb(s):
    # Check for both 'acb' and its reverse forms in the string
    for i in range(len(s)-2):
        substr = s[i:i+3]
        if substr == 'acb' or (substr[::-1] == 'acb' and substr != 'bca'):
            return i
    return -1

def process_string(s):
    steps = [s]
    while True:
        changed = False
        
        # Operation 1: starts with 'bb'
        if s.startswith('bb'):
            s = s[0] + s[2:]
            changed = True
            
        # Operation 2: contains 'acb' or its reverse (not 'bca')
        elif (pos := find_acb(s)) != -1:
            substr = s[pos:pos+3]
            if substr == 'acb':
                s = s[:pos] + 'bca' + s[pos+3:]
            else:  # reverse case
                s = s[:pos] + 'bca' + s[pos+3:]
            changed = True
            
        # Operation 3: ends with 'aa'
        elif s.endswith('aa'):
            s = s[:-2] + 'cc'
            changed = True
            
        # Operation 4: length > 15
        elif len(s) > 15:
            mid = len(s) // 2
            s = s[:mid] + s[mid+1:]
            changed = True
            
        # Operation 5: ends with 'cc'
        elif s.endswith('cc'):
            s = 'a' + s[:-2] + 'b'
            changed = True
            
        # Operation 6: starts with 'aa'
        elif s.startswith('aa'):
            s = s[1:]
            changed = True
            
        if not changed:
            break
            
        steps.append(s)
    
    print("Steps:")
    for i, step in enumerate(steps):
        print(f"Step {i}: {step}")
    return s

initial = "babcbcaccacaaac"
final = process_string(initial)
print(f"\nFinal result: {final}")