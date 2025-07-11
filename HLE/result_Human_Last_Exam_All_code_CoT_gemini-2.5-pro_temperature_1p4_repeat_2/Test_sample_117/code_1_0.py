def find_pskov_voivode_1700():
    """
    This function identifies and prints the name of the voivode who succeeded
    Ivan Ivanovich Golovin in Pskov in 1700.
    """
    # Answer choices provided by the user
    choices = {
        'A': 'Boris Petrovich Sheremetev',
        'B': 'Peter Matveevich Apraksin',
        'C': 'Mikhail Mikhailovich Golitsyn',
        'D': 'Vasily Borisovich Bukhvostov',
        'E': 'Fyodor Yuryevich Romodanovsky'
    }

    # Historical research indicates that Ivan Ivanovich Golovin was the voivode in 1699.
    # At the end of 1699, leading into 1700, Tsar Peter the Great prepared for the
    # Great Northern War against Sweden. He appointed his trusted field marshal,
    # Boris Petrovich Sheremetev, as the commander-in-chief of the main Russian army
    # group, with headquarters in Pskov. As part of this powerful military appointment,
    # Sheremetev also assumed the role of voivode of Pskov in 1700.
    correct_choice_key = 'A'

    successor = choices[correct_choice_key]

    print(f"The successor to Ivan Ivanovich Golovin as the voivode of Pskov in 1700 was:")
    print(f"({correct_choice_key}) {successor}")

find_pskov_voivode_1700()