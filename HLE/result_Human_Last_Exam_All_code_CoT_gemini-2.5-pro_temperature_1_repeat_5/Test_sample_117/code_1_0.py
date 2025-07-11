def solve_historical_question():
    """
    This function identifies and prints the answer to a historical question.
    """
    question = "Who became the next voivode after Ivan Ivanovich Golovin in the city of Pskov in 1700?"
    
    options = {
        "A": "Boris Petrovich Sheremetev",
        "B": "Peter Matveevich Apraksin",
        "C": "Mikhail Mikhailovich Golitsyn",
        "D": "Vasily Borisovich Bukhvostov",
        "E": "Fyodor Yuryevich Romodanovsky"
    }

    # After the defeat at Narva in November 1700, Tsar Peter the Great
    # appointed Boris Petrovich Sheremetev as the commander-in-chief of the
    # Russian army based in Pskov to defend against the Swedes.
    # He effectively took over the military and administrative control of the city.
    correct_answer_key = "A"
    correct_answer_name = options[correct_answer_key]

    print("Question: " + question)
    print("\nAnswer Choices:")
    for key, value in options.items():
        print(f"{key}. {value}")
    
    print("\n-------------------------------------------")
    print(f"The correct answer is: {correct_answer_key}. {correct_answer_name}")
    print("-------------------------------------------")

solve_historical_question()