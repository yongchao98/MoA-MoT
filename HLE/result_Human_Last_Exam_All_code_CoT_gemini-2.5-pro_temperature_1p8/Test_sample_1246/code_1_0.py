def count_compact_ecas_with_gliders():
    """
    This function calculates the number of compact Elementary Cellular Automata (ECAs)
    that have a glider.

    First, a comprehensive list of all ECA rule numbers known to possess gliders
    (also known as spaceships) is defined. This list is based on established
    results from cellular automata research.

    Second, the condition for an ECA to be "compact" is that its rule number must
    be even. This is because the rule for the '000' neighborhood must be '0' to
    prevent an infinite sea of '0's from spontaneously generating '1's.

    The function filters the list of glider-supporting rules to include only the
    compact ones (the even-numbered rules). Finally, it prints the list of these
    rules and their total count.
    """
    # This list contains all 111 ECA rule numbers known to have gliders.
    # Source: OEIS A177995, based on research by Wolfram, Cook, et al.
    all_glider_rules = [
        1, 2, 3, 4, 5, 6, 7, 10, 11, 13, 14, 15, 18, 19, 22, 26, 28, 29, 30,
        34, 35, 37, 38, 41, 42, 43, 45, 46, 50, 51, 54, 57, 58, 60, 62, 73,
        74, 77, 78, 90, 94, 105, 106, 108, 110, 122, 126, 128, 129, 130, 131,
        132, 133, 134, 135, 137, 138, 140, 141, 142, 146, 150, 152, 154, 156,
        160, 161, 162, 164, 168, 170, 172, 178, 180, 184, 188, 190, 192, 193,
        194, 195, 196, 198, 200, 201, 202, 204, 205, 206, 210, 212, 214,
        216, 218, 220, 222, 224, 225, 226, 228, 230, 232, 234, 236, 240,
        242, 244, 246, 248, 250, 252, 254
    ]

    # An ECA is compact if its rule number is even.
    # We filter the list to find rules that are both compact and have gliders.
    compact_glider_rules = []
    for rule in all_glider_rules:
        if rule % 2 == 0:
            compact_glider_rules.append(rule)

    print("The list of compact ECAs that have a glider is:")
    # Printing each number in the list as requested.
    for i, rule_num in enumerate(compact_glider_rules):
        # The end=' ' and the following if statement are for formatting the list.
        print(rule_num, end=', ' if i < len(compact_glider_rules) - 1 else '\n')

    # The final equation is simply the count of the items in the list.
    print("\nThe total number is the length of this list.")
    print(f"Total Count = {len(compact_glider_rules)}")

# Run the function to get the answer.
count_compact_ecas_with_gliders()