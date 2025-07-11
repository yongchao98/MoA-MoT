def find_voivode_successor():
    """
    This script identifies the successor to Ivan Ivanovich Golovin
    as the voivode of Pskov in 1700.
    """

    # The historical question is about the voivode of Pskov in 1700.
    # The predecessor was Ivan Ivanovich Golovin (term ended in 1699).
    # The year 1700 marks the beginning of the Great Northern War, and Pskov
    # was a crucial military headquarters.

    # Historical research indicates Boris Petrovich Sheremetev was appointed
    # as the voivode and commander-in-chief in Pskov for the start of the war.

    answer_choices = {
        'A': 'Boris Petrovich Sheremetev',
        'B': 'Peter Matveevich Apraksin',
        'C': 'Mikhail Mikhailovich Golitsyn',
        'D': 'Vasily Borisovich Bukhvostov',
        'E': 'Fyodor Yuryevich Romodanovsky'
    }

    correct_answer_key = 'A'
    correct_name = answer_choices[correct_answer_key]

    print("Question: Who became the next voivode after Ivan Ivanovich Golovin in the city of Pskov in 1700?")
    print("\nAnswer Choices:")
    for key, value in answer_choices.items():
        print(f"{key}. {value}")

    print("\n----------------------------------------")
    print(f"The correct answer is: {correct_answer_key}. {correct_name}")
    print("Boris Petrovich Sheremetev was a Field Marshal appointed by Peter the Great as commander-in-chief and voivode in Pskov at the start of the Great Northern War in 1700.")

find_voivode_successor()