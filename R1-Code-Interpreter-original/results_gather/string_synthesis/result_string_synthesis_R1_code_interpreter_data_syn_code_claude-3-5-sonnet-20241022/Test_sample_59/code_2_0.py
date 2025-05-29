def analyze_possibilities():
    blocks = {
        '[A]': 0, '[B]': 3, '[C]': 0,
        '{A}': 3, '{B}': 0, '{C}': 0,
        '(A)': 0, '(B)': 0, '(C)': 0
    }
    
    print("Starting state:")
    for block, count in sorted(blocks.items()):
        if count > 0:
            print(f"{count}{block}", end=" ")
    print("\n")
    
    # Try Method 2 ([A] + [B] -> {C}) - Not possible, no [A]
    # Try Method 3 ([B] + [C] -> {B}) - Not possible, no [C]
    # Try Method 4 ([C] + [C] -> {C}) - Not possible, no [C]
    # Try Method 5 ({A} + {C} -> (A) + (B)) - Not possible, no {C}
    # Try Method 6 ({B} + {B} -> (C)) - Not possible, no {B}
    
    print("No further combinations possible")
    print("\nFinal state:")
    for block, count in sorted(blocks.items()):
        if count > 0:
            print(f"{count}{block}", end=" ")

analyze_possibilities()