def find_relevant_sections():
    """
    This function identifies and prints the sections of the UK's Employment Rights Act 1996
    relevant to the procedure for amending section 43FA(3)(c) via statutory instrument.
    """
    # Section to be amended, providing the context for the action.
    target_section = "43FA(3)(c)"

    # Section establishing that the power to make regulations is exercised by statutory instrument.
    instrument_type_section = "236(1)"

    # Section dictating the parliamentary procedure (affirmative resolution) for any statutory
    # instrument that amends the primary Act.
    procedure_section = "236(5)"

    # Combine the sections into a comma-separated list.
    relevant_sections = [target_section, instrument_type_section, procedure_section]
    
    # Print the result as a comma-separated string.
    print(','.join(relevant_sections))

find_relevant_sections()