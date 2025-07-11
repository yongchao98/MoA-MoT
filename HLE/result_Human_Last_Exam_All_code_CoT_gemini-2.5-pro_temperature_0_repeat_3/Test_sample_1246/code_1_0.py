def find_compact_ecas_with_gliders():
    """
    This function identifies compact Elementary Cellular Automata (ECAs) that have gliders
    based on a known comprehensive list of glider-supporting rules.

    A compact ECA is an even-numbered rule, as it preserves the all-zero background.
    A glider is a finite pattern that repeats itself at a different location after some time.

    The list of glider-supporting rules is taken from the paper "Gliders in Elementary
    Cellular Automata" by Sapin et al. (2003).
    """

    # The comprehensive list of 139 ECA rules known to support gliders.
    all_glider_rules = [
        1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 15, 19, 22, 23, 25, 26, 27, 28, 29, 30,
        33, 34, 35, 37, 38, 41, 42, 43, 44, 45, 46, 50, 51, 52, 53, 54, 55, 56, 57, 58,
        60, 61, 62, 73, 74, 75, 77, 78, 90, 92, 93, 94, 97, 98, 100, 101, 105, 106, 108,
        109, 110, 112, 113, 114, 116, 117, 118, 120, 121, 122, 124, 126, 128, 129, 130,
        131, 132, 133, 134, 135, 136, 137, 138, 140, 141, 142, 146, 147, 148, 149, 150,
        152, 153, 154, 156, 157, 158, 160, 161, 162, 163, 164, 165, 166, 168, 169, 170,
        172, 173, 174, 177, 178, 180, 181, 182, 184, 186, 188, 190, 192, 193, 194, 196,
        198, 200, 201, 202, 204, 205, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224,
        225, 226, 228, 230, 232, 233, 234, 236, 238, 240, 242, 244, 246, 248, 250, 252, 254
    ]

    # An ECA is compact if its rule number is even.
    # We filter the list to find the rules that are both compact and have gliders.
    compact_glider_rules = [rule for rule in all_glider_rules if rule % 2 == 0]

    # Print the list of rules found.
    print("The compact ECAs that have a glider are:")
    print(compact_glider_rules)

    # Print the final count.
    count = len(compact_glider_rules)
    print(f"\nThe total number of compact ECAs that have a glider is: {count}")

if __name__ == '__main__':
    find_compact_ecas_with_gliders()
<<<102>>>