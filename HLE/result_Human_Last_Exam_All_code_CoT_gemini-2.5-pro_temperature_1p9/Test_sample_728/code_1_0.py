import collections

# A dictionary to store the films and their Bechdel Test rating.
# A rating of 3 means the film passes all three criteria.
# The structure is {letter: (Title, Year, Bechdel Score)}
film_data = {
    'a': ('Igla', 1988, 3),
    'b': ('Ghost Town', 2008, 3),
    'c': ('Girls Will Be Girls', 2003, 3),
    'd': ('War Dogs', 2016, 1),
    'e': ('Slither', 2006, 3),
    'f': ('John Dies at the End', 2012, 3),
    'g': ('Man Who Knew Too Much, The', 1934, 2),
    'h': ('Ladies In Retirement', 1941, 3),
    'i': ('Internship, The', 2013, 1),
    'j': ('Tinker Bell and the Lost Treasure', 2009, 3)
}

# A list to store the letters of the films that pass the test.
passing_films_letters = []

# Sort the films by their letter (a, b, c...) for consistent output.
sorted_films = collections.OrderedDict(sorted(film_data.items()))

print("Evaluating films based on the Bechdel Test:")
print("A film passes if it meets all three criteria (Score = 3/3).\n")

# Iterate through the films and check their score.
for letter, data in sorted_films.items():
    title, year, score = data
    if score == 3:
        passing_films_letters.append(letter)
        print(f"- Film '{title}' ({year}) passes.")
    else:
        print(f"- Film '{title}' ({year}) fails.")

# Format the final answer as a comma-separated string.
final_answer = ",".join(passing_films_letters)

print("\nThe letters corresponding to the films that pass the Bechdel Test are:")
print(final_answer)
