def solve_sequence():
    """
    Solves the number sequence puzzle by finding the nested patterns.
    """
    # Step 1: Define the initial sequence
    s = [2, 11, 23, 51, 119]
    print(f"The initial sequence is: {', '.join(map(str, s))}")
    print("Step 1: Find the pattern x_n = 2 * x_{n-1} + c_n")
    
    # Step 2: Derive the sequence of constants 'c'
    c = []
    for i in range(len(s) - 1):
        constant = s[i+1] - 2 * s[i]
        c.append(constant)
        print(f"{s[i+1]} = 2 * {s[i]} + {constant}")
    
    print("\nThe sequence of constants 'c' is:", ', '.join(map(str, c)))
    print("\nStep 2: Find the pattern in the constants 'c' by looking at the differences 'd'.")

    # Step 3: Derive the sequence of differences 'd'
    d = []
    for i in range(len(c) - 1):
        diff = c[i+1] - c[i]
        d.append(diff)
        print(f"{c[i+1]} - {c[i]} = {diff}")
        
    print("\nThe sequence of differences 'd' is:", ', '.join(map(str, d)))
    print("\nStep 3: Find the pattern in the differences 'd'.")
    
    # Step 4: Establish the pattern in 'd' and predict the next term
    # d[1] * 3 = d[2] -> 4 * 3 = 12
    ratio = d[2] // d[1]
    print(f"The pattern in 'd' is multiplication by {ratio} (from the 2nd term onwards): {d[1]} * {ratio} = {d[2]}")
    next_d = d[-1] * ratio
    print(f"The next difference is: {d[-1]} * {ratio} = {next_d}")
    
    # Step 5: Predict the next constant 'c'
    next_c = c[-1] + next_d
    print(f"The next constant 'c' is: {c[-1]} + {next_d} = {next_c}")
    
    # Step 6: Calculate the final answer
    next_s = 2 * s[-1] + next_c
    print(f"\nFinally, the next term in the sequence is calculated:")
    print(f"Next Term = 2 * {s[-1]} + {next_c} = {next_s}")

solve_sequence()
<<<A>>>