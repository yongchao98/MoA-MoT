import sys

def find_relevant_sections():
    """
    This function identifies and prints the relevant sections of the UK's Employment Rights Act 1996.

    The logic for identifying the sections is as follows:
    1.  The power to amend Part IVA of the Act, which includes section 43FA, is granted by section 43K. Specifically, s43K(1)(c) allows for amending any provision of that Part.
    2.  Section 236(1) of the Act specifies that any power for the Secretary of State to make an order is "exercisable by statutory instrument".
    3.  Section 236(5) dictates the parliamentary procedure. It mandates the stricter 'affirmative procedure' for orders made under section 43K if they meet the conditions in the subsequent subsection.
    4.  Section 236(6) defines those conditions. An adjustment to s43FA(3)(c) would alter the definition of a protected disclosure, thus meeting the condition in s236(6) and triggering the affirmative procedure required by s236(5).

    These four subsections are therefore all relevant for determining the procedure.
    """
    
    # List of relevant subsections identified through legal analysis
    # s43K(1)(c): The power to amend the relevant part of the Act.
    # s236(1): The requirement to use a statutory instrument.
    # s236(5): The section specifying the affirmative procedure for certain s43K orders.
    # s236(6): The condition which triggers the affirmative procedure.
    relevant_subsections = "43K(1)(c),236(1),236(5),236(6)"

    # Print the comma-separated list of subsections as requested
    print(relevant_subsections)

# Execute the function to find and print the relevant sections.
find_relevant_sections()
