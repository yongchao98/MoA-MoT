# Define the genres of the given movies
movies_genres = {
    "Stargate": ["Science Fiction", "Adventure", "Action"],
    "The Fugitive": ["Thriller", "Action"],
    "Pulp Fiction": ["Crime", "Drama"],
    "Dances with Wolves": ["Historical", "Drama"]
}

# Define the genres of the options
options_genres = {
    "Forrest Gump": ["Drama", "Comedy"],
    "Kevin Smith Too Fat For 40": ["Comedy"],
    "Scorched": ["Comedy"],
    "Stand and Deliver": ["Drama"]
}

# Function to find the most similar movie based on genre overlap
def find_similar_movie(movies_genres, options_genres):
    max_overlap = 0
    most_similar = None
    for option, genres in options_genres.items():
        overlap = sum(1 for genre in genres if any(genre in movie_genres for movie_genres in movies_genres.values()))
        if overlap > max_overlap:
            max_overlap = overlap
            most_similar = option
    return most_similar

# Find the most similar movie
similar_movie = find_similar_movie(movies_genres, options_genres)
print(similar_movie)