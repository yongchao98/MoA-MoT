def find_successor():
    """
    This function identifies and prints the name of the voivode who succeeded
    Ivan Ivanovich Golovin in Pskov in the specified year.
    """
    city = "Pskov"
    year = 1700
    predecessor = "Ivan Ivanovich Golovin"
    
    # Based on historical records, Boris Petrovich Sheremetev was appointed
    # as the voivode of Pskov to command the army group there at the start
    # of the Great Northern War.
    successor = "Boris Petrovich Sheremetev"
    
    print(f"The person who became the next voivode after {predecessor} in the city of {city} in {year} was:")
    print(successor)

find_successor()