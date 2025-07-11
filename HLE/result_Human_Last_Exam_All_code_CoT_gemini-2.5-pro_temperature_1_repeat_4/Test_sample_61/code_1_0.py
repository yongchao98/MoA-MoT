def find_relevant_sections():
    """
    This function identifies and prints the sections of the Employment Rights Act 1996
    that determine the parliamentary procedure for regulations made under section 43FA(3)(c).

    - 43FA(3)(c): Grants the specific power to the Secretary of State to make regulations.
    - 236(1): Confirms that any power to make regulations under the Act is exercisable by statutory instrument.
    - 236(2): Establishes the default procedure for such statutory instruments as the negative resolution procedure ("subject to annulment").
    - 236(5): Lists the specific regulation-making powers that are subject to the affirmative procedure. This section is relevant because the power under 43FA(3)(c) is NOT listed, confirming the default negative procedure from 236(2) applies.
    """
    
    # The relevant sections, when read together, determine the procedure.
    relevant_sections = [
        "43FA(3)(c)",
        "236(1)",
        "236(2)",
        "236(5)"
    ]
    
    # Format the output as a comma-separated list
    output = ",".join(relevant_sections)
    
    print(output)

find_relevant_sections()