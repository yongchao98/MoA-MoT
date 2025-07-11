def find_relevant_sections():
    """
    This function identifies and prints the sections of the UK's Employment Rights Act 1996
    that determine the parliamentary procedure for statutory instruments made to adjust
    section 43FA(3)(c).

    The process is as follows:
    1.  Section 43FA(3)(c): This subsection grants the power to the Secretary of State
        to make regulations specifying a "person of a description". This is the starting point.
    2.  Section 236: This is the general section in the Act that governs how all orders and
        regulations are made.
    3.  Section 236(1): This subsection mandates that any power for the Secretary of State
        to make regulations under the Act must be exercised "by statutory instrument".
    4.  Section 236(3): This subsection specifies the parliamentary procedure. It states that
        any statutory instrument for which the affirmative procedure is not explicitly required
        "shall be subject to annulment in pursuance of a resolution of either House of
        Parliament". This is the definition of the negative procedure. Since the power under
        43FA(3)(c) is not listed as requiring the affirmative procedure elsewhere in section 236,
        this default rule applies.

    Therefore, the complete chain of relevant provisions is 43FA(3)(c) (the power),
    236(1) (the instrument type), and 236(3) (the procedure).
    """
    
    # The section granting the power to make regulations.
    power_section = "43FA(3)(c)"
    
    # The section specifying that the power is exercised by statutory instrument.
    instrument_type_section = "236(1)"
    
    # The section specifying the negative procedure for this type of instrument.
    procedure_section = "236(3)"
    
    # Combine the relevant sections into a comma-separated list for output.
    # The prompt asks to "output each number in the final equation", which is interpreted
    # as printing each section identifier in the final list.
    final_list = f"{power_section},{instrument_type_section},{procedure_section}"
    
    print(final_list)

find_relevant_sections()