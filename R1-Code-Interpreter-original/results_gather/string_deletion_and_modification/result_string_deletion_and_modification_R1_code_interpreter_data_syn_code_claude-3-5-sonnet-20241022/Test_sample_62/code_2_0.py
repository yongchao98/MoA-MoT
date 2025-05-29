def apply_rules(s):
    print(f"Starting with: {s}")
    
    # Count b's and check if rule 2 applies
    b_count = s.count('b')
    print(f"Number of b's: {b_count}")
    
    if b_count > 0 and b_count % 2 == 0:
        result = s + 'ab'
        print(f"Rule 2 applies: {result}")
    else:
        result = s
        print("Rule 2 doesn't apply")
        
    return result

# Test the string
s = "bbbcbbaaabbbb"
final = apply_rules(s)
print(f"\nFinal result: {final}")