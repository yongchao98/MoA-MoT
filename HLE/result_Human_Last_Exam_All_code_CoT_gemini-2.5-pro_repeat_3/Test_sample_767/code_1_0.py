def solve_limit_problem():
    """
    This function calculates the value of the limit lim_{N->inf} p(N)/N.

    As explained in the plan, the value of this limit is equal to the number of
    non-negative integers n for which the Fibonacci number F_n is less than or
    equal to 25. This is because these are the only solutions that contribute
    a term proportional to N to the total count of solutions p(N).
    """

    # The maximum value F_n can take, derived from the constraint on g.
    limit_f = 25

    print(f"We are looking for the number of integers n >= 0 such that F_n <= {limit_f}.")
    print("This number will be the value of the limit.")
    print("-" * 20)

    count = 0
    n = 0
    f_prev = 0  # Represents F_0
    f_curr = 1  # Represents F_1

    # We handle n=0 separately as the loop starts from F_1.
    f_n = f_prev
    if f_n <= limit_f:
        print(f"n = {n}: F_{n} = {f_n}. This is a valid solution (since {f_n} <= {limit_f}).")
        count += 1
    n += 1

    # Now we loop through Fibonacci numbers starting from F_1 until they exceed the limit.
    f_n = f_curr
    while f_n <= limit_f:
        print(f"n = {n}: F_{n} = {f_n}. This is a valid solution (since {f_n} <= {limit_f}).")
        count += 1
        
        # Calculate the next Fibonacci number
        f_next = f_prev + f_curr
        f_prev = f_curr
        f_curr = f_next
        f_n = f_curr # The new f_n for the next iteration check
        n += 1

    print(f"n = {n}: F_{n} = {f_n}. This is not a valid solution (since {f_n} > {limit_f}).")
    print("-" * 20)
    
    final_answer = count
    print(f"The total number of valid n is {final_answer}.")
    print(f"Therefore, the final equation for the limit is:")
    print(f"lim_{{N->inf}} p(N)/N = {final_answer}")

solve_limit_problem()