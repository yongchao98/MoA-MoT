def transform_string(s):
    seen_states = set()
    while s not in seen_states:
        seen_states.add(s)
        original_s = s  # Keep track of the original state before applying any rule
        print(f"Current string: {s}")  # Debug: print the current state
        
        if s.startswith('cc'):
            s = s[2:]
            print("Applied Rule 1: Remove 'cc' from start")
        elif s.startswith('bb'):
            s = s[0] + s[2:]
            print("Applied Rule 2: Remove second 'b' from start")
        elif s.endswith('ac'):
            s = s[:-2] + 'cb'
            print("Applied Rule 3: Replace 'ac' at end with 'cb'")
        elif s.startswith('cb'):
            s = 'aa' + s[2:-1]
            print("Applied Rule 4: Replace 'cb' at start with 'aa' and remove last character")
        elif s.startswith('ac'):
            s = 'zz' + s[2:]
            print("Applied Rule 5: Replace 'ac' at start with 'zz'")
        
        # If no rule was applied, break the loop
        if s == original_s:
            print("No rule applied, breaking the loop")
            break

    return s

initial_string = "aaabacbcaccc"
final_string = transform_string(initial_string)
print(f"Final string: {final_string}")