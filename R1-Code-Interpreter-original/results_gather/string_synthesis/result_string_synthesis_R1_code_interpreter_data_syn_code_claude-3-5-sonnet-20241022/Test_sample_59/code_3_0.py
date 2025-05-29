def verify_final_state():
    # Initialize current state
    state = {
        '[A]': 0, '[B]': 3, '[C]': 0,
        '{A}': 3, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    print("Checking all possible rule applications:")
    
    # Check each rule's requirements and current state
    rules = [
        ("[A]+[B]+[C]->{A}", lambda s: s['[A]']>=1 and s['[B]']>=1 and s['[C]']>=1),
        ("[A]+[B]->{C}", lambda s: s['[A]']>=1 and s['[B]']>=1),
        ("[B]+[C]->{B}", lambda s: s['[B]']>=1 and s['[C]']>=1),
        ("[C]+[C]->{C}", lambda s: s['[C]']>=2),
        ("{A}+{C}->(A)+(B)", lambda s: s['{A}']>=1 and s['{C}']>=1),
        ("{B}+{B}->(C)", lambda s: s['{B}']>=2)
    ]
    
    for rule_name, check in rules:
        if check(state):
            print(f"Rule {rule_name} can be applied")
        else:
            print(f"Rule {rule_name} cannot be applied - missing required blocks")
    
    print("\nCurrent state:")
    for block, count in sorted(state.items()):
        if count > 0:
            print(f"{count}{block}", end=" ")

verify_final_state()