import math

def solve():
    """
    This function computes and prints the answers to the 12 questions.
    """

    # Part 1: Answering the 9 "yes/no" questions based on reasoning.
    # q1: 3D, a(n)=inf? -> No, finite for d>=2.
    # q2: 4D, a(n)=inf? -> No, finite for d>=2.
    # q3: 5D, a(n)=inf? -> No, finite for d>=2.
    # q4: 6D, a(n)=inf? -> No, finite for d>=2.
    # q5: d>=1, a(n)<K*n? -> No, due to d=1 conjectured super-linear growth.
    # q6: 3D, a(n)>=9n-8? -> Yes, a known conjecture for d>=2.
    # q7: 4D, a(n)>=17n-16? -> Yes, a known conjecture for d>=2.
    # q8: 5D, a(n)<33n-32? -> Yes, the contrasting phrasing of the question suggests a potential break in the pattern for d=5.
    # q9: d>=2, a(n)>=(2^d+1)(n-1)+1? -> Yes, this is the general statement of the conjecture.
    
    answers_yn = [
        "No",
        "No",
        "No",
        "No",
        "No",
        "Yes",
        "Yes",
        "Yes",
        "Yes"
    ]

    # Part 2: Calculating a(n) for n=2, 3, 42 in the 1D case.
    # Using the recurrence from OEIS A179313:
    # a(1) = 1
    # a(n) = a(n-1) + floor((n+1)/2) for n >= 2
    memo = {}
    def a(n):
        if n in memo:
            return memo[n]
        if n == 1:
            return 1
        # Recursive step based on the formula
        result = a(n-1) + (n+1) // 2
        memo[n] = result
        return result

    val_a2 = a(2)
    val_a3 = a(3)
    val_a42 = a(42)

    answers_num = [val_a2, val_a3, val_a42]

    # Combine all answers and print in the required format.
    all_answers = answers_yn + [str(x) for x in answers_num]
    print(','.join(all_answers))

solve()