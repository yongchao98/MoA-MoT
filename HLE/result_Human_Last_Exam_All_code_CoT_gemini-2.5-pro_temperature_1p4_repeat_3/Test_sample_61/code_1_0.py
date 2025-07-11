def find_relevant_sections():
    """
    Identifies and prints the sections of the Employment Rights Act 1996 relevant 
    to the procedure for making regulations under section 43FA(3)(c).
    """
    
    # Section 43FA(3)(c): Grants the power to the Secretary of State to make regulations.
    section_power = "43FA(3)(c)"
    
    # Section 236(1): Specifies that this power must be exercised by statutory instrument.
    section_instrument_type = "236(1)"
    
    # Section 236(2)(b): Specifies the parliamentary procedure (negative resolution)
    # for such a statutory instrument.
    section_procedure = "236(2)(b)"

    # Combine the sections into a comma-separated list as requested.
    relevant_sections_list = [section_power, section_instrument_type, section_procedure]
    
    print(','.join(relevant_sections_list))

find_relevant_sections()