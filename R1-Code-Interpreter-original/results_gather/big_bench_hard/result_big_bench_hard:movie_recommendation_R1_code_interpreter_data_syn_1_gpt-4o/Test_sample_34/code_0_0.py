# Define the genres of the given movies
given_movies_genres = {
    "The Sixth Sense": ["supernatural", "thriller", "mystery", "drama"],
    "The Lord of the Rings: The Two Towers": ["fantasy", "adventure"],
    "The Usual Suspects": ["crime", "thriller", "mystery"],
    "Forrest Gump": ["drama", "romance", "historical"]
}

# Define the genres of the options
options_genres = {
    "Nothing to Declare": ["comedy"],
    "Crossfire": ["film noir", "crime", "drama"],
    "End of Days": ["supernatural", "action", "horror"],
    "Star Wars Episode IV - A New Hope": ["science fiction", "fantasy", "adventure"]
}

# Function to count genre matches
def count_genre_matches(given_genres, option_genres):
    matches = 0
    for genre in option_genres:
        for movie_genres in given_genres.values():
            if genre in movie_genres:
                matches += 1
    return matches

# Calculate matches for each option
matches = {option: count_genre_matches(given_movies_genres, genres) for option, genres in options_genres.items()}

# Find the option with the maximum matches
most_similar_movie = max(matches, key=matches.get)
print(most_similar_movie)