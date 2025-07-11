def find_relevant_sections():
    """
    This function identifies and prints the sections of the UK's Employment Rights Act 1996
    relevant for determining the statutory instrument procedure for adjusting section 43FA(3)(c).
    """
    # Section 43FA(3)(c): Grants the power to the Secretary of State to make an order.
    # Section 236(1): Specifies that this power is exercised by statutory instrument.
    # Section 236(2): Sets the default 'negative procedure' for such instruments.
    # Section 236(5): Lists the exceptions requiring the 'affirmative procedure'. An order
    # under s.43FA is not listed, confirming the default procedure applies.
    relevant_sections = [
        "43FA(3)(c)",
        "236(1)",
        "236(2)",
        "236(5)"
    ]
    
    print(",".join(relevant_sections))

find_relevant_sections()