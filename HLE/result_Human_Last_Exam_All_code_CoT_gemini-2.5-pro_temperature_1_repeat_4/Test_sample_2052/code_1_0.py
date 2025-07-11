def solve_question():
    """
    This function provides the answers to the user's questions based on the theory of matrix models.

    (a) Superintegrability for n=3: The matrix model for n=3 is equivalent to the Kontsevich model,
        which is a well-known superintegrable system. So the answer is Yes.

    (b) Necessity of W_{-4} for n=4: For n>=4, the partition function defined by the character expansion
        is a "non-standard" solution to the Virasoro constraints, which is not equivalent to the
        result from the matrix integral. The operator W_{-4} is introduced in the literature
        specifically to generate this non-standard solution. Therefore, it is considered necessary.
        So the answer is Yes.
    """
    answer_a = "Yes"
    answer_b = "Yes"
    
    print(f"(a) {answer_a}; (b) {answer_b}")

solve_question()