def find_relevant_sections():
    """
    Identifies and prints the relevant sections of the Employment Rights Act 1996
    for determining the statutory instrument procedure for adjusting section 43FA(3)(c).
    """
    # The target section to be amended.
    target_section = "43FA(3)(c)"

    # Section 236(2) specifies that the power is exercisable by statutory instrument.
    instrument_procedure_section = "236(2)"

    # Section 236(3) specifies the default parliamentary procedure (negative procedure),
    # as a power to amend 43FA is not listed in 236(5) as requiring the affirmative procedure.
    parliamentary_procedure_section = "236(3)"

    # Combine the sections into a comma-separated list.
    relevant_sections_list = [
        target_section,
        instrument_procedure_section,
        parliamentary_procedure_section
    ]

    final_output = ",".join(relevant_sections_list)
    print("The relevant sections of the Employment Rights Act 1996 are:")
    print(final_output)

find_relevant_sections()