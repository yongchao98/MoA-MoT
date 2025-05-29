# Define the genres of the given movies
given_movies_genres = {
    "Star Wars Episode VII - The Force Awakens": ["Science Fiction", "Adventure", "Action"],
    "Inside Out": ["Animation", "Family", "Adventure", "Fantasy"],
    "Pulp Fiction": ["Crime", "Drama", "Dark Comedy"],
    "Raiders of the Lost Ark": ["Action", "Adventure", "Fantasy"]
}

# Define the genres of the options
options_genres = {
    "Ernest Rides Again": ["Comedy", "Adventure"],
    "Forrest Gump": ["Drama", "Romance", "Comedy"],
    "The Proposal": ["Romantic Comedy"],
    "Everything You Always Wanted to Know About Sex But Were Afraid to Ask": ["Comedy"]
}

# Function to find the best match based on common genres
def find_best_match(given_movies_genres, options_genres):
    common_genres = set()
    for genres in given_movies_genres.values():
        common_genres.update(genres)
    
    best_match = None
    max_common = 0
    
    for movie, genres in options_genres.items():
        common_count = len(set(genres) & common_genres)
        if common_count > max_common:
            max_common = common_count
            best_match = movie
    
    return best_match

# Find the best match
best_match = find_best_match(given_movies_genres, options_genres)
print(best_match)