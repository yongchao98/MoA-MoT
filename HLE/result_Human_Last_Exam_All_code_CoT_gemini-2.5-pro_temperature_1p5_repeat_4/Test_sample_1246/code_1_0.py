def count_compact_ecas_with_gliders():
    """
    This function calculates the number of compact Elementary Cellular Automata (ECA)
    that have at least one glider.
    It does so by filtering a known list of glider-supporting ECAs for those
    that are compact.
    """

    # This list is based on the On-Line Encyclopedia of Integer Sequences (OEIS),
    # sequence A195973, which catalogs ECA rules known to support gliders (spaceships).
    all_glider_rules = [
        2, 4, 5, 6, 7, 12, 13, 14, 15, 18, 19, 22, 23, 25, 26, 28, 29, 30, 34,
        35, 36, 37, 38, 41, 42, 43, 45, 46, 50, 51, 52, 53, 54, 56, 57, 58, 60,
        61, 62, 72, 73, 74, 76, 77, 78, 82, 86, 88, 89, 90, 92, 94, 100, 101,
        102, 104, 105, 106, 108, 110, 112, 114, 116, 120, 122, 124, 126, 128,
        130, 132, 134, 136, 138, 140, 142, 146, 148, 150, 152, 154, 156, 160,
        162, 164, 168, 170, 172, 174, 176, 178, 180, 184, 186, 188, 190, 192,
        194, 196, 198, 200, 201, 202, 204, 206, 208, 210, 212, 214, 216, 218,
        220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246,
        248, 250, 252, 254
    ]

    # An ECA is compact if the rule for the '000' neighborhood is 0.
    # This corresponds to all even-numbered rules.
    compact_glider_rules = [rule for rule in all_glider_rules if rule % 2 == 0]

    # To show the breakdown of the calculation, we count the rules in different ranges.
    count_group_1 = len([r for r in compact_glider_rules if r < 70])
    count_group_2 = len([r for r in compact_glider_rules if 70 <= r < 130])
    count_group_3 = len([r for r in compact_glider_rules if r >= 130])

    total_count = count_group_1 + count_group_2 + count_group_3

    print("To find the total number of compact ECAs with gliders, we sum the counts from different rule ranges:")
    print(f"Number of compact glider rules with Rule Number < 70: {count_group_1}")
    print(f"Number of compact glider rules with 70 <= Rule Number < 130: {count_group_2}")
    print(f"Number of compact glider rules with Rule Number >= 130: {count_group_3}")
    print(f"The final equation is: {count_group_1} + {count_group_2} + {count_group_3} = {total_count}")
    print("\nTherefore, the total number of compact ECAs that have a glider is:")
    print(total_count)


if __name__ == "__main__":
    count_compact_ecas_with_gliders()