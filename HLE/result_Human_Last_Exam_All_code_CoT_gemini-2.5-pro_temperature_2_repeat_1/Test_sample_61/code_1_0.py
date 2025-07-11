import sys

def find_relevant_sections():
    """
    This function identifies and prints the relevant sections of the UK's Employment Rights Act 1996.

    - Section 43FA(3)(c): Grants the power to the Secretary of State to specify, by order made by statutory instrument, matters which constitute a "relevant breach".
    - Section 236(3): Sets out the default parliamentary procedure for statutory instruments made under the Act. It states they are "subject to annulment in pursuance of a resolution of either House of Parliament" (the negative procedure), unless specified otherwise. Section 43FA is not listed as an exception.
    - Section 236(5): Lists the powers for which an order requires a draft to be "approved by a resolution of each House of Parliament" (the affirmative procedure). The power in section 43FA is not included in this list, which confirms by exclusion that the procedure in section 236(3) applies.

    Therefore, these three sections are key to determining the correct procedure.
    """
    section_granting_power = "43FA(3)(c)"
    section_defining_negative_procedure = "236(3)"
    section_defining_affirmative_procedure = "236(5)"

    # The final list of relevant sections, comma-separated
    final_list = f"{section_granting_power},{section_defining_negative_procedure},{section_defining_affirmative_procedure}"
    
    print(final_list)

find_relevant_sections()