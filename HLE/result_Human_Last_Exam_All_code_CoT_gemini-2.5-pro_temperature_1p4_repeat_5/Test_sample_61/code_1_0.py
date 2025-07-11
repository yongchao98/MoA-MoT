def find_relevant_sections():
    """
    Identifies and prints the sections of the Employment Rights Act 1996
    relevant for determining the procedure for an order under section 43FA(3)(c).
    """
    # Section 43FA(3)(c) contains the power to make the order.
    # Section 236(1) specifies this order must be a statutory instrument.
    # Section 236(2) specifies the parliamentary procedure (negative procedure)
    # as 43FA is not listed as an exception.
    relevant_sections = [
        "43FA(3)(c)",
        "236(1)",
        "236(2)"
    ]

    # Print the sections as a comma-separated list
    print(','.join(relevant_sections))

find_relevant_sections()