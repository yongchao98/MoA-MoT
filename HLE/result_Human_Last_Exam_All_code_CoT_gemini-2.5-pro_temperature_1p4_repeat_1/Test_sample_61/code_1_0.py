def find_relevant_sections():
    """
    Identifies and prints the sections of the Employment Rights Act 1996
    relevant for determining the procedure for amending section 43FA(3)(c).
    """
    # Section 43FA(5): Grants the power to amend section 43FA(3) by order.
    # Section 43FA(6): Specifies the exact procedure (affirmative resolution) for that order.
    # Section 236(3): Defines the affirmative procedure, providing context.
    # Section 236(5): Defines the default negative procedure, providing contrast and context.
    relevant_sections = [
        "43FA(5)",
        "43FA(6)",
        "236(3)",
        "236(5)"
    ]
    print(",".join(relevant_sections))

find_relevant_sections()