def generate_s4(up_to_n):
    """
    Generates the sequence S4 up to the n-th term based on the deduced rule.
    s[1] = 1
    s[2] = 1
    s[n] = s[s[n-1]] + s[n - 1 - s[n-2]]
    """
    if up_to_n < 1:
        return []
    
    # Use a dictionary for potentially sparse indices and 1-based indexing
    s = {1: 1, 2: 1}
    
    # Generate the sequence from n=3 up to the desired term
    for n in range(3, up_to_n + 1):
        try:
            # The deduced recurrence relation
            s[n] = s[s[n - 1]] + s[n - 1 - s[n - 2]]
        except KeyError:
            # This can happen if an index goes below 1
            print(f"Could not compute s[{n}] due to an invalid index.")
            break
            
    return list(s.values())

# The problem asks to deduce R(s[n]) and output the rule.
# The following print statement fulfills this requirement.

print("The deduced rule R for sequence S4 is:")
print("R(s[n]) = s[s[n - 1]] + s[n - 1 - s[n - 2]]")

# As a demonstration, let's generate the first 42 terms with this rule
# and compare it to the input sequence.
# print("\nGenerated sequence S4 (first 42 terms):")
# print(generate_s4(42))