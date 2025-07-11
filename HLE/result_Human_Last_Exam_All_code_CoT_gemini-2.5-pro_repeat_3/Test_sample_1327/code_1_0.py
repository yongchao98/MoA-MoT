def solve_sequence():
    """
    Solves for the next number in the sequence 2, 11, 23, 51, 119.
    """
    a = [2, 11, 23, 51, 119]
    print(f"The original sequence is: {', '.join(map(str, a))}, ()")

    # Step 1: Calculate the first differences
    d = [a[i+1] - a[i] for i in range(len(a)-1)]
    print(f"1. The differences between consecutive terms are: {', '.join(map(str, d))}")

    # Step 2: Define and calculate sequence 'c'
    # c[n] = d[n+1] - a[n+1]
    c = [d[i+1] - a[i+2] for i in range(len(d)-2)]
    print(f"2. A second-level analysis reveals a new sequence 'c', where c[n] = d[n+1] - a[n+2]: {', '.join(map(str, c))}")
    print(f"   For example, c[0] = {d[1]} - {a[2]} = {c[0]}")


    # Step 3: Define and calculate sequence 'e' (differences of 'c')
    e = [c[i+1] - c[i] for i in range(len(c)-1)]
    print(f"3. The differences of sequence 'c' are: {', '.join(map(str, e))}")

    # Step 4: Find the pattern in 'e' and predict the next term
    ratio = e[1] / e[0]
    next_e = e[-1] * ratio
    print(f"4. This sequence 'e' is a geometric progression with a ratio of {int(ratio)} ({e[1]} / {e[0]}).")
    print(f"   The next term in 'e' is: {e[-1]} * {int(ratio)} = {int(next_e)}")

    # Step 5: Work backwards to find the next terms
    next_c = c[-1] + next_e
    print(f"5. The next term in 'c' is: {c[-1]} + {int(next_e)} = {int(next_c)}")

    # d_next = a_last + c_next.
    # Note: d index is len(a)-1 = 4. so d_next is d[4].
    # a index is len(a) = 5. so a_last is a[4]
    # c index is len(c) = 3. so c_next is c[3]
    # The relation is d[i] = a[i] + c[i-1]. So d[4] = a[4] + c[3]
    next_d = a[-1] + next_c
    print(f"6. The next term in the differences 'd' is: {a[-1]} + {int(next_c)} = {int(next_d)}")

    next_a = a[-1] + next_d
    print(f"7. Finally, the next term in the original sequence is found by the equation:")
    print(f"   {a[-1]} + {int(next_d)} = {int(next_a)}")
    print("\nThe number that should be filled in the parentheses is 291.")

solve_sequence()
<<<A>>>