def solve_sequence():
    """
    Solves the number sequence puzzle and explains the steps.
    """
    # The given sequence
    a = [2, 11, 23, 51, 119]
    print(f"The given sequence is: {', '.join(map(str, a))}, ...")

    print("\nStep 1: Find the relationship between terms. We test a pattern of the form a_n+1 = 2 * a_n + c_n, starting from n=2.")
    
    # Calculate the sequence c_n
    c2 = a[2] - 2 * a[1]
    print(f"For n=2: {a[2]} = 2 * {a[1]} + c_2  =>  c_2 = {a[2]} - {2 * a[1]} = {c2}")
    
    c3 = a[3] - 2 * a[2]
    print(f"For n=3: {a[3]} = 2 * {a[2]} + c_3  =>  c_3 = {a[3]} - {2 * a[2]} = {c3}")
    
    c4 = a[4] - 2 * a[3]
    print(f"For n=4: {a[4]} = 2 * {a[3]} + c_4  =>  c_4 = {a[4]} - {2 * a[3]} = {c4}")
    
    c_sequence = [c2, c3, c4]
    print(f"\nStep 2: We have a new sequence c = {', '.join(map(str, c_sequence))}. Let's find the pattern here by looking at the differences.")
    
    # Calculate the differences of the c_sequence
    diff1 = c_sequence[1] - c_sequence[0]
    print(f"First difference: {c_sequence[1]} - {c_sequence[0]} = {diff1}")
    
    diff2 = c_sequence[2] - c_sequence[1]
    print(f"Second difference: {c_sequence[2]} - {c_sequence[1]} = {diff2}")
    
    print("\nThe differences (4, 12) form a geometric progression with a common ratio of 3 (12 / 4 = 3).")

    # Predict the next difference and the next term in c_sequence
    next_diff = diff2 * 3
    print(f"\nStep 3: Predict the next difference: {diff2} * 3 = {next_diff}")
    
    c5 = c_sequence[-1] + next_diff
    print(f"Predict the next term in sequence c: {c_sequence[-1]} + {next_diff} = {c5}")
    
    # Calculate the final answer
    next_a = 2 * a[-1] + c5
    print(f"\nStep 4: Calculate the missing number in the original sequence using a_6 = 2 * a_5 + c_5.")
    print(f"The final equation is: 2 * {a[-1]} + {c5} = {next_a}")
    
solve_sequence()