def solve_and_print_answers():
    """
    Solves a series of 12 theoretical questions about a grid growth model
    and prints the answers as a comma-separated string.
    """

    # --- Reasoning for Higher Dimensions (d>=2) ---
    # A key result for this problem states that for any dimension d>=2,
    # if we start with n>=2 ones, it is possible to choose an initial
    # configuration and a growth strategy (e.g., a spiral or helical path)
    # such that the process never terminates. This means a_d(n) = infinity for d>=2, n>=2.
    # For n=1, we can only place '1's, so a_d(1) = 1 in any dimension.

    answers = []

    # Q1: In 3d, is it true that a(n)=infinity if n is large?
    # Yes. "n is large" means we can take n>=2.
    answers.append("Yes")

    # Q2: The same question in 4d?
    # Yes, for the same reason.
    answers.append("Yes")

    # Q3: The same question in 5d?
    # Yes, for the same reason.
    answers.append("Yes")

    # Q4: The same question in 6d?
    # Yes, for the same reason.
    answers.append("Yes")

    # Q5: In d>=1 dimension is it true that a(n)<K*n?
    # No. For d>=2 and n>=2, a(n) is infinite and cannot be bounded by a finite value K*n.
    answers.append("No")

    # The next questions introduce a potential bound B_d(n) = (2^d+1)(n-1)+1.
    # Q6: In 3d, is it true that a(n) >= 9n-8 for every n? (Note 9 = 2^3+1)
    # Yes. For n>=2, a(n) is infinite, which is greater than any finite number.
    # For n=1, a(1)=1, and the inequality is 1 >= 9(1)-8=1, which is true.
    answers.append("Yes")

    # Q7: On the 4d case, is it true that a(n) >= 17n-16 for every n? (Note 17 = 2^4+1)
    # Yes, for the same reason as in 3D.
    answers.append("Yes")

    # Q8: On the 5d case, is it true that a(n) < 33n-32 if n is large? (Note 33 = 2^5+1)
    # No. For large n (n>=2), a(n) is infinite and cannot be less than a finite number.
    answers.append("No")

    # Q9: In general for d>=2 dimension is it true that a(n) >= (2^d+1)(n-1)+1?
    # Yes. This generalizes the logic from Q6 and Q7. Since a_d(n) is infinite for
    # n>=2, the inequality holds. For n=1, it simplifies to 1>=1.
    answers.append("Yes")


    # --- Reasoning for 1D Case ---
    # In 1D, a cell p has only two neighbors: p-1 and p+1.
    # To place a number k>1 at p, the cells p-1 and p+1 must already be occupied.
    # This implies that p itself is "shielded" and no other empty cell can have p as a neighbor.
    # Therefore, any number k>1 placed on the grid can never contribute to future sums.
    # All subsequent numbers must be formed by summing the initial '1's. The only possible sum is 1+1=2.
    # So, we can place 2 (if n>=2), but then we can't place 3. The process stops.
    # Thus, a(n) = 1 for n=1, and a(n) = 2 for n>=2.

    # Q10: On the 1d case what is a(2)?
    # For n=2, we can place 2, so a(2)=2.
    answers.append("2")

    # Q11: On the 1d case what is a(3)?
    # For n=3, a(3)=2.
    answers.append("2")

    # Q12: On the 1d case what is a(42)?
    # For n=42, a(42)=2.
    answers.append("2")

    print(",".join(answers))

solve_and_print_answers()