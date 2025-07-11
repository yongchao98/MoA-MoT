def solve():
    """
    Calculates the sum of squares of coefficients for the given polynomial expansion.
    Let S(n) be the sum of the squares of the coefficients of P_n(x) = prod_{i=0 to n} (1+x^(3^i)+x^(2*3^i)+x^(3*3^i)).
    The sequence S(n) follows the recurrence relation:
    S(n) = 10*S(n-1) - 26*S(n-2) + 8*S(n-3) for n >= 3.
    """
    # Initial values for the recurrence S(0), S(1), S(2).
    s = [4, 22, 124]

    print("The recurrence relation is S(n) = 10*S(n-1) - 26*S(n-2) + 8*S(n-3)")
    print("Initial values:")
    print(f"S(0) = {s[0]}")
    print(f"S(1) = {s[1]}")
    print(f"S(2) = {s[2]}")
    print("-" * 20)
    print("Calculating S(n) up to n=19:")

    # We need to find S(19).
    # The loop calculates S(n) for n from 3 to 19.
    for n in range(3, 20):
        s_prev_1 = s[n-1]
        s_prev_2 = s[n-2]
        s_prev_3 = s[n-3]
        
        s_n = 10 * s_prev_1 - 26 * s_prev_2 + 8 * s_prev_3
        s.append(s_n)
        
        # Outputting each number in the equation for each step
        print(f"S({n}) = 10 * {s_prev_1} - 26 * {s_prev_2} + 8 * {s_prev_3} = {s_n}")

    final_answer = s[19]
    print("-" * 20)
    print(f"The final result for the sum of squares of the coefficients, S(19), is: {final_answer}")

solve()