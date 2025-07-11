import sys

def find_relevant_sections():
    """
    This function outlines the reasoning and identifies the relevant sections of the
    UK Employment Rights Act 1996 for determining the statutory instrument procedure
    for adjusting section 43FA(3)(c).
    """

    # Step 1: The target for adjustment is section 43FA(3)(c), which relates to
    # the definition of a 'worker' for the purposes of whistleblowing protection
    # in Part IVA of the Act.

    # Step 2: The most likely power to make such a change is found in section 43K.
    # Section 43K(1) gives the Secretary of State the power to make an order
    # to modify the application of this Part of the Act.
    power_section = "43K(1)"

    # Step 3: The general rules for making statutory instruments under the Act are
    # in section 236 ("Orders and regulations").
    # Section 236(1) states that these powers must be exercised by statutory instrument.
    mechanism_section = "236(1)"

    # Step 4: Section 236(5) specifies the precise parliamentary procedure for
    # certain orders. It explicitly states that a statutory instrument containing an
    # order under section 43K "shall not be made unless a draft of the instrument
    # has been laid before Parliament and approved by a resolution of each House
    # of Parliament." This is the affirmative procedure.
    procedure_section = "236(5)"

    # Step 5: The final list consists of the section granting the specific power
    # and the sections defining the mechanism and procedure for that power.
    relevant_sections = [power_section, mechanism_section, procedure_section]

    # Print the comma-separated list as requested.
    print(",".join(relevant_sections))

# Execute the function to find and print the sections.
find_relevant_sections()

# <<<43K(1),236(1),236(5)>>>