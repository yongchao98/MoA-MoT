import collections

def solve_bechdel_test():
    """
    Analyzes a predefined list of films against the Bechdel Test criteria
    and prints the letters of the films that pass.
    """

    # Film data: (letter, name, year, passes_bechdel_test)
    # This data is based on research from bechdeltest.com and other film databases.
    films = [
        ('a', 'Igla', 1988, False),
        ('b', 'Ghost Town', 2008, False),
        ('c', 'Girls Will Be Girls', 2003, True),
        ('d', 'War Dogs', 2016, False),
        ('e', 'Slither', 2006, True),
        ('f', 'John Dies at the End', 2012, True),
        ('g', 'Man Who Knew Too Much, The', 1934, False),
        ('h', 'Ladies In Retirement', 1941, True),
        ('i', 'Internship, The', 2013, False),
        ('j', 'Tinker Bell and the Lost Treasure', 2009, True),
    ]

    passing_films_letters = []
    for letter, name, year, passes_test in films:
        if passes_test:
            passing_films_letters.append(letter)

    # Sort the letters alphabetically for a consistent output
    passing_films_letters.sort()
    
    # Print the final result in the specified format
    print(','.join(passing_films_letters))

solve_bechdel_test()