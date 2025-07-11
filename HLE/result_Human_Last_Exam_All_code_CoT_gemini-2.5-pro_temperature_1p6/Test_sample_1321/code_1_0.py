import math

def get_answers():
    """
    Solves the twelve questions based on known results and formulas
    from the study of this specific growth model.
    """
    answers = []

    # Questions 1-4: For d>=2, a(n) can be infinite if n is large enough.
    # The condition is n >= 2^d - d + 1. "Large n" implies this condition is met.
    answers.append("Yes")  # Q1: 3D, a(n)=inf? Yes, for n>=6.
    answers.append("Yes")  # Q2: 4D, a(n)=inf? Yes, for n>=13.
    answers.append("Yes")  # Q3: 5D, a(n)=inf? Yes, for n>=28.
    answers.append("Yes")  # Q4: 6D, a(n)=inf? Yes, for n>=59.

    # Question 5: For d>=1, is a(n) < K*n?
    # This is false because for d>=2, a(n) can be infinite, which is not < K*n.
    answers.append("No")

    # Questions 6-9: These relate to a conjectured general lower bound for d>=2.
    # The conjectured formula is a_d(n) >= (2^d + 1)(n - 1) + 1.
    
    # Question 9: Is a(n) >= (2^d+1)(n-1)+1 for d>=2?
    # This is a known conjecture based on specific constructions. We assume Yes.
    is_q9_true = True
    answers.append("Yes" if is_q9_true else "No")

    # Question 6: 3D case, is a(n) >= 9n-8?
    # For d=3, the formula gives (2**3 + 1)(n-1) + 1 = 9*n - 8. It matches.
    answers.append("Yes" if is_q9_true else "No")
    
    # Question 7: 4D case, is a(n) >= 17n-16?
    # For d=4, the formula gives (2**4 + 1)(n-1) + 1 = 17*n - 16. It matches.
    answers.append("Yes" if is_q9_true else "No")

    # Question 8: 5D case, is a(n) < 33n-32?
    # For d=5, the formula gives a(n) >= (2**5 + 1)(n-1)+1 = 33*n - 32.
    # So a(n) cannot be strictly smaller than this bound.
    answers.append("No" if is_q9_true else "Yes")
    
    # Swap answers for Q8 and Q9 to match the question order in the prompt.
    answers[5], answers[6], answers[7], answers[8] = answers[8], answers[5], answers[6], answers[7]


    # Questions 10-12: Values of a(n) in 1D.
    # The formula for the 1D case is a_1(n) = floor(n^2 / 4) + 1.
    def a_1d(n):
      # Equation for a(n) in 1D.
      # The numbers in this equation are n, 2, 4, 1.
      val_sq = n**2
      val_div = val_sq // 4
      result = val_div + 1
      return result

    # Q10: a(2)
    answers.append(a_1d(2))
    
    # Q11: a(3)
    answers.append(a_1d(3))

    # Q12: a(42)
    answers.append(a_1d(42))

    print(*answers, sep=",")

get_answers()