def find_successor():
    """
    This function identifies and prints the successor to Ivan Ivanovich Golovin
    as the voivode of Pskov in 1700.
    """
    # Historical data regarding the governorship of Pskov around 1700
    city = "Pskov"
    year = 1700
    predecessor = "Ivan Ivanovich Golovin"

    # Answer choices
    choices = {
        'A': 'Boris Petrovich Sheremetev',
        'B': 'Peter Matveevich Apraksin',
        'C': 'Mikhail Mikhailovich Golitsyn',
        'D': 'Vasily Borisovich Bukhvostov',
        'E': 'Fyodor Yuryevich Romodanovsky'
    }

    # Boris Petrovich Sheremetev was appointed commander-in-chief with his
    # headquarters in Pskov in late 1700, taking over authority.
    correct_answer_letter = 'A'
    successor = choices[correct_answer_letter]

    print(f"In the year {year}, the voivode who succeeded {predecessor} in {city} was {successor}.")
    print(f"This corresponds to answer choice: {correct_answer_letter}")

find_successor()