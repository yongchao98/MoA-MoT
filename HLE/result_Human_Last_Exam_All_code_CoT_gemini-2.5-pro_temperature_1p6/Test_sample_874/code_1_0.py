def solve():
    """
    Finds the specific tuple (a, b, c, d) and computes the required expression.
    """
    limit = 10_000_000

    # Step 1 & 2: Generate Tribonacci numbers
    T = [0, 0, 1]
    while T[-1] <= limit:
        next_t = T[-1] + T[-2] + T[-3]
        T.append(next_t)

    # Step 3: Find the largest n for the construction
    best_n = 0
    for n in range(1, len(T) // 3):
        if 3 * n + 1 >= len(T):
            break
        
        # Following the construction formulas
        t_3n = T[3 * n]
        t_3n_minus_1 = T[3 * n - 1]
        t_3n_minus_2 = T[3 * n - 2]
        t_3n_minus_3 = T[3 * n - 3]
        t_3n_plus_1 = T[3*n+1]
        
        # d_n = t_{3n-1}
        d = t_3n_minus_1
        # c_n = T_{3n} - T_{3n-3}
        c = t_3n - t_3n_minus_3
        # b_n = t_{3n}
        b = t_3n
        # a_n = T_{3n+1} - T_{3n-1}
        a = t_3n_plus_1 - t_3n_minus_1
        
        if a <= limit and b <= limit and c <= limit and d <= limit:
            best_n = n
        else:
            break

    # Step 4: Construct the final tuple using the best n
    n = best_n
    t_3n = T[3 * n]
    t_3n_minus_1 = T[3 * n - 1]
    t_3n_minus_2 = T[3 * n - 2]
    t_3n_minus_3 = T[3 * n - 3]
    t_3n_plus_1 = T[3*n+1]

    d = t_3n_minus_1
    c = t_3n - t_3n_minus_3
    b = t_3n
    a = t_3n_plus_1 - t_3n_minus_1

    print(f"The largest n that fits the constraints is n = {n}.")
    print(f"The tuple (a,b,c,d) found is:")
    print(f"a = {t_3n_plus_1} - {t_3n_minus_1} = {a}")
    print(f"b = {t_3n} = {b}")
    print(f"c = {t_3n} - {t_3n_minus_3} = {c}")
    print(f"d = {t_3n_minus_1} = {d}")
    print(f"The sum a+b+c+d is {a+b+c+d}.")

    # Step 5: Compute the expression
    result = (a + b - c - d) % 1000
    
    print("\nCalculating (a + b - c - d) mod 1000:")
    print(f"= ({a} + {b} - {c} - {d}) mod 1000")
    print(f"= (({a % 1000}) + ({b % 1000}) - ({c % 1000}) - ({d % 1000})) mod 1000")
    print(f"= ({a % 1000} + {b % 1000} - {c % 1000} - {d % 1000}) mod 1000")
    a_mod = a % 1000
    b_mod = b % 1000
    c_mod = c % 1000
    d_mod = d % 1000
    
    print(f"= {a_mod} + {b_mod} - {c_mod} - {d_mod} (mod 1000)")
    sum_ab = (a_mod + b_mod)
    sub_cd = (c_mod + d_mod)
    print(f"= {sum_ab} - {sub_cd} (mod 1000)")
    
    final_val = (sum_ab - sub_cd) % 1000
    print(f"= {sum_ab - sub_cd} (mod 1000)")
    print(f"= {final_val}")
    
    return final_val

final_answer = solve()
print(f"\nThe final computed expression a + b - c - d mod 1000 is {final_answer}.")
print(f"<<<{final_answer}>>>")
