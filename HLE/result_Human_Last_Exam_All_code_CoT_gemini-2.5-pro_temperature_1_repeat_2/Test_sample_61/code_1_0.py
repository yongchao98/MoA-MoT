def find_relevant_sections():
    """
    This function identifies and prints the sections of the UK's Employment Rights Act 1996
    that are relevant for determining the procedure for adjusting section 43FA(3)(c)
    via a statutory instrument.

    - Section 43KA(1): This subsection grants the Secretary of State the specific power by order to amend section 43FA(3)(c). This is the foundational enabling power.
    - Section 43KA(3): This subsection specifies that any order made under this section must be made by statutory instrument.
    - Section 43KA(4): This subsection defines the specific parliamentary procedure required, stating that a draft of the instrument must be laid before and approved by a resolution of each House of Parliament (the affirmative resolution procedure).
    """
    
    # The relevant section that grants the power to amend section 43FA(3)(c)
    power_granting_section = "43KA(1)"
    
    # The subsection specifying the instrument type
    instrument_type_subsection = "43KA(3)"
    
    # The subsection specifying the parliamentary procedure
    procedure_subsection = "43KA(4)"
    
    # Combine the relevant subsections into a comma-separated list
    relevant_sections_list = f"{power_granting_section},{instrument_type_subsection},{procedure_subsection}"
    
    # Print the final comma-separated list
    print(relevant_sections_list)

find_relevant_sections()