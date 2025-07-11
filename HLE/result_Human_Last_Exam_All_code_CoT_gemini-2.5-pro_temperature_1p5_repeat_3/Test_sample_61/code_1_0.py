def find_relevant_sections():
    """
    This script identifies and prints the sections of the Employment Rights Act 1996
    that determine the procedure for making a statutory instrument to adjust
    section 43FA(3)(c) of the Act.
    """

    # Step 1: The power to make an order relating to prescribed persons for disclosures
    # is granted by section 43FA(1). This is the specific power being exercised.
    power_section = "43FA(1)"

    # Step 2: Section 236(1) of the Act specifies that any power for the Secretary of
    # State to make an order under the Act is exercisable by statutory instrument.
    instrument_section = "236(1)"

    # Step 3: Section 236(3) establishes the default parliamentary procedure. It states that
    # a statutory instrument under the Act is subject to annulment by a resolution of
    # either House of Parliament (the "negative procedure"). An order under section 43FA
    # is not listed among the exceptions that require the stricter "affirmative procedure".
    procedure_section = "236(3)"

    # Combine the relevant sections into the final list.
    final_list = f"{power_section},{instrument_section},{procedure_section}"

    print(f"The power to create the order is in {power_section}.")
    print(f"The requirement to use a statutory instrument is in {instrument_section}.")
    print(f"The parliamentary procedure (negative resolution) is specified in {procedure_section}.")
    print("\nThe final comma-separated list of relevant sections is:")
    print(final_list)

find_relevant_sections()