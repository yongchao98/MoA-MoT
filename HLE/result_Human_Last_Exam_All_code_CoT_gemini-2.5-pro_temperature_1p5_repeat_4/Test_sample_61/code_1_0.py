def find_relevant_sections():
    """
    Identifies and prints the sections of the UK's Employment Rights Act 1996
    relevant to the procedure for amending section 43FA(3)(c).
    """
    # Section 43FA(3) grants the power to the Secretary of State to make an "order".
    power_section = "43FA(3)"

    # Section 236(1) confirms that an "order" is made by statutory instrument.
    instrument_type_section = "236(1)"

    # Section 236(5) specifies the parliamentary procedure for this statutory instrument
    # (negative resolution), as 43FA is not listed in 236(4) for the affirmative procedure.
    procedure_section = "236(5)"

    relevant_sections = [power_section, instrument_type_section, procedure_section]

    # Print the sections as a comma-separated list.
    print(",".join(relevant_sections))

find_relevant_sections()