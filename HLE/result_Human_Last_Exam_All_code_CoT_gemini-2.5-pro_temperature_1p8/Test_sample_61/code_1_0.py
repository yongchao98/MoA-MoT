def find_relevant_sections():
    """
    Identifies and prints the sections of the UK's Employment Rights Act 1996
    relevant to the procedure for adjusting section 43FA(3)(c).
    """

    # Section 43FA(3)(c): This is the subsection to be adjusted. It contains the
    # specific power for the Secretary of State to make regulations specifying a
    # "description" of person for national security disclosures.
    power_source = "43FA(3)(c)"

    # Section 236(1): This general provision states that any power for the
    # Secretary of State to make regulations under the Act must be exercised
    # "by statutory instrument".
    instrument_type = "236(1)"

    # Section 236(3): This sets the default parliamentary procedure for most
    # statutory instruments under the Act. It specifies the negative resolution
    # procedure ("subject to annulment"), but also notes there are exceptions.
    default_procedure = "236(3)"

    # Section 236(5): This is the crucial exception. It overrides the default
    # procedure for specific regulations, explicitly including those "under section 43FA".
    # It mandates the affirmative resolution procedure, stating that a draft instrument
    # "shall not be made unless a draft of the instrument has been laid before and
    # approved by a resolution of each House of Parliament."
    specific_procedure = "236(5)"

    relevant_sections = [
        power_source,
        instrument_type,
        default_procedure,
        specific_procedure
    ]

    print(','.join(relevant_sections))

find_relevant_sections()