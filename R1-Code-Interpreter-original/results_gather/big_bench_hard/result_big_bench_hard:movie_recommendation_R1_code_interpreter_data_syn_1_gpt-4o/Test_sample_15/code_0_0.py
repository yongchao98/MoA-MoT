# Define the genres and themes of the given movies
given_movies = {
    "Forrest Gump": ["drama", "romance", "comedy"],
    "The Silence of the Lambs": ["thriller", "horror", "crime"],
    "Seven": ["thriller", "crime"],
    "Fargo": ["crime", "dark comedy"]
}

# Define the genres and themes of the options
options = {
    "Gandhi": ["historical", "drama"],
    "Schindler's List": ["historical", "drama"],
    "Dogfight": ["romance", "drama"],
    "Repo Man": ["science fiction", "crime", "dark comedy"]
}

# Function to find the most similar movie
def find_similar_movie(given_movies, options):
    # Count the number of matching genres/themes
    similarity_scores = {}
    for option, genres in options.items():
        score = 0
        for movie, movie_genres in given_movies.items():
            score += len(set(genres) & set(movie_genres))
        similarity_scores[option] = score
    # Find the option with the highest score
    most_similar = max(similarity_scores, key=similarity_scores.get)
    return most_similar

# Find and print the most similar movie
most_similar_movie = find_similar_movie(given_movies, options)
print(most_similar_movie)