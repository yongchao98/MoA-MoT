def find_relevant_sections():
    """
    Identifies and prints the sections of the UK's Employment Rights Act 1996
    relevant for determining the procedure for adjusting section 43FA(3)(c).
    """
    # The section to be adjusted is the starting point.
    target_section = "43FA(3)(c)"

    # Section 236(1) states that powers to make orders are exercisable by statutory instrument.
    power_by_si = "236(1)"

    # Section 236(3) specifies the default "negative procedure" for statutory instruments under the Act.
    negative_procedure = "236(3)"

    # Section 236(5) specifies which orders are subject to the "affirmative procedure".
    # One must check this list to see if the specific power is included.
    affirmative_procedure = "236(5)"

    # These sections are all relevant for determining the correct procedure.
    relevant_sections = [
        target_section,
        power_by_si,
        negative_procedure,
        affirmative_procedure
    ]

    print(",".join(relevant_sections))

find_relevant_sections()