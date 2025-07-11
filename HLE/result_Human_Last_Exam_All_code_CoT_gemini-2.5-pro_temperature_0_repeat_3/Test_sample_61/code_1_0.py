def find_relevant_sections():
    """
    This function identifies and prints the sections of the UK's Employment Rights Act 1996
    relevant for determining the procedure for adjusting section 43FA(3)(c).
    """

    # Section 43FA(4) grants the specific power to the Secretary of State to make regulations
    # for the purposes of section 43FA(3)(c).
    power_section = "43FA(4)"

    # Section 236(1) specifies that any power to make regulations under the Act is
    # exercisable by statutory instrument.
    instrument_type_section = "236(1)"

    # Section 236(3) sets the default parliamentary procedure for statutory instruments
    # under the Act. It states they are "subject to annulment in pursuance of a resolution
    # of either House of Parliament" (the negative procedure), unless they are one of the
    # specified exceptions. The power under 43FA(4) is not listed as an exception.
    procedure_section = "236(3)"

    # The relevant sections are the one granting the power and the ones defining the procedure.
    relevant_sections = [power_section, instrument_type_section, procedure_section]

    # Print the result as a comma-separated list.
    print(','.join(relevant_sections))

find_relevant_sections()