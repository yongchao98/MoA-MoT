def solve_literary_puzzle():
    """
    This script solves the puzzle by explaining the literary pun
    at its core.
    """
    # 1. The puzzle refers to the Vienna Ringstrasse, a ring of boulevards.
    # The key word is "ring".
    
    # 2. In Joseph Brodsky's original Russian essay, he makes a pun.
    # The Russian word for "ring" is 'кольцо'.
    russian_word_for_ring = "кольцо"
    pronunciation_ring = "kol'tso"

    # 3. The surname of the English poet, when transliterated into Russian,
    # starts with the same root as the word for "ring".
    poet_surname_in_russian = "Кольридж"
    pronunciation_poet = "Kol'ridzh"
    
    # 4. The "equation" is the shared root between the two words.
    shared_root = "Коль-"

    print("The puzzle is based on a pun in a Russian essay by Joseph Brodsky.")
    print("He connects the 'ring' of boulevards in Vienna to a poet's name.")
    print("-" * 30)
    print("Here is the logical 'equation':")
    print(f"1. The Russian word for 'ring' is: {russian_word_for_ring} (pronounced '{pronunciation_ring}')")
    print(f"2. The Russian spelling of the poet's name is: {poet_surname_in_russian} (pronounced '{pronunciation_poet}')")
    print(f"3. The connection is the shared root: '{shared_root}'")
    print("-" * 30)

    # 5. This reveals the English surname.
    english_surname = "Coleridge"
    print(f"Therefore, the English poet's surname is: {english_surname}")

solve_literary_puzzle()