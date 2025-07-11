def find_successor():
    """
    This function identifies and prints the name of the person who became the voivode of Pskov after Ivan Ivanovich Golovin in 1700.
    """
    choices = {
        'A': "Boris Petrovich Sheremetev",
        'B': "Peter Matveevich Apraksin",
        'C': "Mikhail Mikhailovich Golitsyn",
        'D': "Vasily Borisovich Bukhvostov",
        'E': "Fyodor Yuryevich Romodanovsky"
    }

    # Historical records indicate Boris Petrovich Sheremetev was appointed commander
    # of the forces in Pskov at the start of the Great Northern War in 1700,
    # making him the de facto successor.
    correct_choice = 'A'
    
    print(f"The person who became the next voivode (military governor) in Pskov in 1700 was: {choices[correct_choice]}")

find_successor()