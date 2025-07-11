def find_relevant_sections():
    """
    This function identifies and prints the relevant sections of the Employment Rights Act 1996.
    
    - 43FA(3)(c): Grants the power to the Secretary of State to make regulations by statutory instrument.
    - 236(1): Confirms that powers to make regulations under the Act are exercised by statutory instrument.
    - 236(2): Specifies the default 'negative procedure' for statutory instruments, which applies to regulations made under 43FA(3)(c) as it is not listed as an exception requiring the affirmative procedure.
    """
    
    sections = [
        "43FA(3)(c)",
        "236(1)",
        "236(2)"
    ]
    
    # The output format is a comma-separated list.
    print(",".join(sections))

find_relevant_sections()