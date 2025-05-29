# List of given movies and their genres
given_movies = {
    "Star Wars Episode V - The Empire Strikes Back": ["science fiction", "action", "adventure"],
    "Independence Day": ["science fiction", "action", "adventure"],
    "Terminator 2: Judgment Day": ["science fiction", "action", "adventure"],
    "Star Wars Episode VI - Return of the Jedi": ["science fiction", "action", "adventure"]
}

# Options and their assumed genres
options = {
    "High School": ["comedy", "drama"],
    "Pushing Tin": ["comedy", "drama"],
    "What Happened": ["unknown"],
    "Miss Simone": ["documentary"],
    "Back to the Future": ["science fiction", "adventure", "comedy"]
}

# Function to find the most similar movie based on genre
def find_similar_movie(given_movies, options):
    # Common genres in given movies
    common_genres = set.intersection(*[set(genres) for genres in given_movies.values()])
    
    # Find the option with the most matching genres
    best_match = None
    max_matches = 0
    for movie, genres in options.items():
        matches = len(common_genres.intersection(genres))
        if matches > max_matches:
            best_match = movie
            max_matches = matches
    
    return best_match

# Find and print the most similar movie
similar_movie = find_similar_movie(given_movies, options)
print(similar_movie)