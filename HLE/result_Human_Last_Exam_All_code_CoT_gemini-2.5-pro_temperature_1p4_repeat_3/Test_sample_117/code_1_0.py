def solve_history_question():
    """
    This function identifies the successor to Ivan Ivanovich Golovin as the voivode of Pskov in 1700.
    """
    question = "Who became the next voivode after Ivan Ivanovich Golovin in the city of Pskov in 1700?"

    choices = {
        'A': 'Boris Petrovich Sheremetev',
        'B': 'Peter Matveevich Apraksin',
        'C': 'Mikhail Mikhailovich Golitsyn',
        'D': 'Vasily Borisovich Bukhvostov',
        'E': 'Fyodor Yuryevich Romodanovsky'
    }

    # Based on historical records, Boris Petrovich Sheremetev was appointed
    # voivode of Pskov and commander of the forces there after the Battle of Narva in late 1700.
    correct_answer_key = 'A'
    correct_answer_name = choices[correct_answer_key]

    print(f"Question: {question}")
    print("Answer Choices:")
    for key, value in choices.items():
        print(f"{key}. {value}")
    
    print("\n----------------------------------------\n")
    print(f"The correct answer is: {correct_answer_key}. {correct_answer_name}")

solve_history_question()