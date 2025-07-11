def find_next_voivode():
    """
    This function solves a historical question by providing the known successor
    to Ivan Ivanovich Golovin as the voivode of Pskov in 1700.
    """
    question = "Who became the next voivode after Ivan Ivanovich Golovin in the city of Pskov in 1700?"

    options = {
        'A': 'Boris Petrovich Sheremetev',
        'B': 'Peter Matveevich Apraksin',
        'C': 'Mikhail Mikhailovich Golitsyn',
        'D': 'Vasily Borisovich Bukhvostov',
        'E': 'Fyodor Yuryevich Romodanovsky'
    }

    # Historical context:
    # In preparation for the Great Northern War, which began in 1700, Tsar Peter the Great
    # appointed Boris Petrovich Sheremetev as the commander-in-chief and voivode
    # of the Russian forces concentrated in Pskov. He arrived in Pskov and took command,
    # effectively succeeding the previous voivode, Ivan Ivanovich Golovin.

    correct_answer_key = 'A'
    correct_answer_name = options[correct_answer_key]

    print(f"Question: {question}")
    print("\nHistorical Analysis:")
    print("In 1700, at the start of the Great Northern War, the role of voivode in a key border city like Pskov was essentially a military governorship.")
    print("Tsar Peter the Great appointed Field Marshal Boris Petrovich Sheremetev to command the main Russian army group based in Pskov.")
    print("Therefore, Boris Petrovich Sheremetev succeeded Ivan Ivanovich Golovin in this role.")
    print(f"\nThe correct answer is: {correct_answer_name}")

find_next_voivode()