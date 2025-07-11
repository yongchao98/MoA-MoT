def find_voivode_successor():
    """
    This function identifies and prints the successor to Ivan Ivanovich Golovin
    as the voivode of Pskov in 1700 from the given options.
    """
    # Historical context: With the start of the Great Northern War in 1700,
    # Tsar Peter I appointed Boris Petrovich Sheremetev as the commander-in-chief
    # and voivode overseeing the region, including Pskov.
    
    # Answer choices provided by the user
    options = {
        'A': 'Boris Petrovich Sheremetev',
        'B': 'Peter Matveevich Apraksin',
        'C': 'Mikhail Mikhailovich Golitsyn',
        'D': 'Vasily Borisovich Bukhvostov',
        'E': 'Fyodor Yuryevich Romodanovsky'
    }
    
    # The correct answer based on historical records
    correct_answer_key = 'A'
    year = 1700
    
    # Print the full answer details
    print(f"The person who became the next voivode in Pskov in the year {year} was:")
    print(f"{correct_answer_key}. {options[correct_answer_key]}")

find_voivode_successor()