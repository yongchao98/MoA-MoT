def find_relevant_sections():
    """
    Identifies and prints the sections of the Employment Rights Act 1996
    relevant to the parliamentary procedure for amending section 43FA(3)(c).
    """
    # Section 43FA(3)(c): Grants the power to the Secretary of State to make an "order".
    section_granting_power = "43FA(3)(c)"

    # Section 236(1): Specifies that an "order" under the Act is a statutory instrument.
    section_defining_instrument = "236(1)"

    # Section 236(3): Specifies the default procedure (negative resolution) for
    # statutory instruments not explicitly listed for the affirmative procedure.
    section_defining_procedure = "236(3)"

    relevant_sections = [
        section_granting_power,
        section_defining_instrument,
        section_defining_procedure
    ]

    print(",".join(relevant_sections))

find_relevant_sections()