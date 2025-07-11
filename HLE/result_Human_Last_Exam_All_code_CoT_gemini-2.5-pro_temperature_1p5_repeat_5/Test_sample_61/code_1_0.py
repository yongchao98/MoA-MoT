def find_relevant_sections():
    """
    Identifies the relevant sections of the Employment Rights Act 1996 for a specific statutory instrument procedure.

    The user is asking which sections determine the procedure for a statutory instrument to adjust section 43FA(3)(c).

    1.  The power to make the statutory instrument (an "order") is granted by section 43FA(3)(c) itself.
        - 43FA(3)(c): "(a relevant member of Parliament is a member of...) any other description of person specified in an order made by the Secretary of State."

    2.  The general rules for how orders under this Act are made are found in section 236.
        - Section 236(1) states that any power to make an order is "exercisable by statutory instrument."

    3.  Section 236(3) sets the default parliamentary procedure for these statutory instruments.
        - It states they are "subject to annulment in pursuance of a resolution of either House of Parliament" (the negative procedure), unless they are one of the specific exceptions listed in 236(3) or 236(5).

    4.  An order under 43FA(3)(c) is not listed as an exception, so the default rule in 236(3) applies.

    Therefore, the relevant sections are 43FA(3)(c), 236(1), and 236(3).
    """
    sections = [
        "43FA(3)(c)",
        "236(1)",
        "236(3)"
    ]
    # Output the result as a comma-separated list
    print(",".join(sections))

find_relevant_sections()