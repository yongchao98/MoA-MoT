def find_relevant_sections():
    """
    Identifies the sections of the Employment Rights Act 1996 relevant for determining
    the parliamentary procedure for statutory instruments.
    """
    # Section 236(3) establishes that powers are exercised by statutory instrument.
    # Section 236(4) establishes the default negative procedure.
    # Section 236(5) lists the exceptions requiring the affirmative procedure.
    # These sections together provide the framework for determining the procedure.
    relevant_sections = ["236(3)", "236(4)", "236(5)"]
    print(','.join(relevant_sections))

find_relevant_sections()