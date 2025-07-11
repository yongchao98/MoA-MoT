def solve_insect_sexing():
    """
    Solves the insect sex identification puzzle based on morphological traits.

    The options are:
    1) M, M
    2) F, F
    3) M, F
    4) F, M
    """

    # Pair A: Cuckoo Bees (e.g., Coelioxys)
    # The left bee has a sharply pointed abdomen, which is a female ovipositor.
    # The right bee has a blunter abdomen with prongs, characteristic of a male.
    # Conclusion: Female, Male (F, M)
    answer_A = 4

    # Pair B: Paper Wasps (Polistes)
    # The left wasp has long, curled antennae, which is a male characteristic.
    # The right wasp has shorter, straight antennae, characteristic of a female.
    # Conclusion: Male, Female (M, F)
    answer_B = 3

    # Pair C: Long-horned Bees (e.g., Melissodes)
    # The left bee has very long antennae and a yellow face patch (clypeus), which are male traits.
    # The right bee has short antennae and a dark face, characteristic of a female.
    # Conclusion: Male, Female (M, F)
    answer_C = 3

    # Format the final answer as requested
    final_answer = f"{answer_A},{answer_B},{answer_C}"
    print(final_answer)

solve_insect_sexing()