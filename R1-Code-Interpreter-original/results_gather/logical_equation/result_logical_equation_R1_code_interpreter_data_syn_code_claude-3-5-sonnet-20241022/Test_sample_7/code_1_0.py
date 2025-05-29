def verify_and_print_details(solution):
    A, B, C, D, E, F, G, H, I = solution
    
    # Print all relationships for verification
    print(f"Solution being tested: {solution}")
    print(f"Checking D({D}) + H({H}) = 13: {D + H == 13}")
    print(f"Checking E({E}) - C({C}) = 15: {E - C == 15}")
    print(f"Checking C({C}) + H({H}) = 5: {C + H == 5}")
    print(f"Checking H({H}) = 4.0C({C}): {H == 4.0 * C}")
    print(f"Checking H({H}) - B({B}) = -117: {H - B == -117}")
    print(f"Checking F({F}) = 3.0C({C}): {F == 3.0 * C}")
    print(f"Checking F({F}) + H({H}) = 7: {F + H == 7}")
    print(f"Checking all numbers are unique: {len(set(solution)) == 9}")
    print(f"Checking all numbers are in valid set: {all(x in [1,3,4,9,16,27,36,80,121] for x in solution)}")
    print("---")

# Test all solutions from the output
solutions = [
    [27, 121, 1, 9, 16, 3, 36, 4, 80],
    [27, 121, 1, 9, 16, 3, 80, 4, 36],
    [36, 121, 1, 9, 16, 3, 27, 4, 80],
    [36, 121, 1, 9, 16, 3, 80, 4, 27],
    [80, 121, 1, 9, 16, 3, 27, 4, 36],
    [80, 121, 1, 9, 16, 3, 36, 4, 27]
]

for solution in solutions:
    verify_and_print_details(solution)