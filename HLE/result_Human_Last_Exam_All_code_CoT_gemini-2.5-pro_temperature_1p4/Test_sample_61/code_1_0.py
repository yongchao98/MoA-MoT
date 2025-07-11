def find_relevant_sections():
    """
    This function identifies and prints the sections of the Employment Rights Act 1996
    that are relevant for determining the procedure for a statutory instrument
    adjusting section 43FA(3)(c).
    
    - 43FA(3)(c): The substantive provision creating the need for a statutory instrument.
    - 236(1): Confirms that powers under the Act are exercised by statutory instrument.
    - 236(2): Establishes the default negative resolution procedure.
    - 236(3), 236(4), 236(5): List the exceptions requiring the affirmative procedure. They must be
      checked to confirm the default procedure applies.
    """
    
    # The substantive provision
    source_provision = "43FA(3)(c)"
    
    # The general procedural provisions in section 236
    procedure_si = "236(1)"
    procedure_negative = "236(2)"
    procedure_affirmative_exceptions_1 = "236(3)"
    procedure_affirmative_exceptions_2 = "236(4)"
    procedure_affirmative_exceptions_3 = "236(5)"
    
    relevant_sections = [
        source_provision,
        procedure_si,
        procedure_negative,
        procedure_affirmative_exceptions_1,
        procedure_affirmative_exceptions_2,
        procedure_affirmative_exceptions_3
    ]
    
    print(",".join(relevant_sections))

find_relevant_sections()