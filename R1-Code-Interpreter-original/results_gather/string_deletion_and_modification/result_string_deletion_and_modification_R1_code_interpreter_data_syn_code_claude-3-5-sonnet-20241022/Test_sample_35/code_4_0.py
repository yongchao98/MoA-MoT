def print_state(s, operation=""):
    print(f"{operation}Current string: {s}")
    print(f"Length: {len(s)}")
    print(f"First char: {s[0]}")
    print(f"Last chars: {s[-2:]}")
    print(f"b count: {s.count('b')}")
    if operation:
        print("-" * 40)

def transform_string(s):
    original = s
    print(f"Starting with: {s}\n")
    
    while True:
        changed = False
        
        # Operation 1: If contains 'ca' (not at start)
        if 'ca' in s[1:]:
            pos = s[1:].find('ca') + 1  # add 1 because we searched from index 1
            s_new = s[:pos] + s[pos+2:]
            print(f"Operation 1: Found 'ca' at position {pos}, removing...")
            print_state(s_new, "After Op1: ")
            s = s_new
            changed = True
            continue

        # Operation 2: Even number of b's
        b_count = s.count('b')
        if b_count > 0 and b_count % 2 == 0:
            s_new = s + 'ab'
            print(f"Operation 2: Found {b_count} b's (even), appending 'ab'")
            print_state(s_new, "After Op2: ")
            s = s_new
            changed = True
            continue

        # Operation 3: Contains 'bca'
        if 'bca' in s:
            pos = s.find('bca')
            s_new = s[:pos] + s[pos+3:]
            print(f"Operation 3: Found 'bca' at position {pos}, removing...")
            print_state(s_new, "After Op3: ")
            s = s_new
            changed = True
            continue

        # Operation 4: Prefix 'cb'
        if s.startswith('cb'):
            s_new = 'aa' + s[2:-1]
            print("Operation 4: Found 'cb' prefix, replacing...")
            print_state(s_new, "After Op4: ")
            s = s_new
            changed = True
            continue

        # Operation 5: Prefix 'ab'
        if s.startswith('ab'):
            s_new = 'ca' + s[2:]
            print("Operation 5: Found 'ab' prefix, replacing...")
            print_state(s_new, "After Op5: ")
            s = s_new
            changed = True
            continue

        # Operation 6: Suffix 'cc'
        if s.endswith('cc'):
            s_new = 'a' + s[:-2] + 'b'
            print("Operation 6: Found 'cc' suffix, replacing...")
            print_state(s_new, "After Op6: ")
            s = s_new
            changed = True
            continue

        if not changed:
            break

    print(f"\nFinal string: {s}")
    print(f"No more operations possible.")
    return s

s = "acaccbcbaacccbcccaac"
transform_string(s)