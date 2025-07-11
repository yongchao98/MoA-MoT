def find_relevant_sections():
    """
    This function identifies and prints the sections of the UK's Employment Rights Act 1996
    that are relevant for determining the parliamentary procedure for a statutory instrument
    adjusting section 43FA(3)(c).
    """
    
    # The section to be adjusted, which defines the subject matter.
    target_section = "43FA(3)(c)"
    
    # The section listing powers subject to the affirmative procedure.
    affirmative_procedure_section = "236(2)"
    
    # The section providing the default negative procedure.
    negative_procedure_section = "236(3)"
    
    # The section listing other powers subject to the affirmative procedure.
    other_affirmative_procedure_section = "236(5)"
    
    # Combine the sections into a comma-separated list for the output.
    relevant_sections_list = [
        target_section,
        affirmative_procedure_section,
        negative_procedure_section,
        other_affirmative_procedure_section
    ]
    
    # Print the final comma-separated list.
    print(','.join(relevant_sections_list))

find_relevant_sections()