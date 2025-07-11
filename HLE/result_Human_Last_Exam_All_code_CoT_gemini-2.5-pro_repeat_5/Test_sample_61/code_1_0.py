def find_relevant_sections():
    """
    Identifies and prints the sections of the UK's Employment Rights Act 1996
    relevant for determining the procedure for adjusting section 43FA(3)(c).
    """
    # The power to adjust s.43FA(3)(c) is in s.43FA(4).
    # The general rules for the procedure are in s.236.
    # s.236(1) confirms it's a statutory instrument.
    # s.236(3) sets the parliamentary procedure (negative procedure) as it's the default
    # and an order under 43FA(4) is not listed as an exception.
    relevant_sections = "43FA(4),236(1),236(3)"
    print(relevant_sections)

find_relevant_sections()