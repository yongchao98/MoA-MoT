def solve_grid_questions():
    """
    This function formulates the answers to the twelve questions and prints them.
    """
    answers = [
        # Q1: 3d, a(n)=inf for large n?
        "Yes",
        # Q2: 4d, a(n)=inf for large n?
        "Yes",
        # Q3: 5d, a(n)=inf for large n?
        "Yes",
        # Q4: 6d, a(n)=inf for large n?
        "Yes",
        # Q5: d>=1, a(n)<K*n?
        "No",
        # Q6: 3d, a(n)>=9n-8?
        "Yes",
        # Q7: 4d, a(n)>=17n-16?
        "Yes",
        # Q8: 5d, a(n)<33n-32?
        "No",
        # Q9: d>=2, a(n)>=(2^d+1)(n-1)+1?
        "Yes",
        # Q10: 1d, a(2)?
        2,
        # Q11: 1d, a(3)?
        4,
        # Q12: 1d, a(42)?
        # For n>=3 in 1D, a(n) follows the formula 3n - 5.
        # For n = 42:
        # 3 * 42 - 5 = 121
        121
    ]

    # The final answer should be a comma-separated string.
    print(",".join(map(str, answers)))

solve_grid_questions()