# Define a dictionary with movie titles and their genres
movies_genres = {
    "The Shawshank Redemption": ["Drama"],
    "Forrest Gump": ["Drama", "Comedy"],
    "Dances with Wolves": ["Drama", "Western"],
    "Mr Holland's Opus": ["Drama", "Music"],
    "PK": ["Comedy", "Drama"],
    "The Great Gatsby": ["Drama", "Romance"],
    "Philadelphia": ["Drama"],
    "Casper Meets Wendy": ["Family", "Fantasy"]
}

# Define the target genres based on the given movies
target_genres = ["Drama"]

# Function to find the most similar movie based on genre
def find_similar_movie(options, target_genres):
    for option in options:
        genres = movies_genres.get(option, [])
        if any(genre in genres for genre in target_genres):
            return option
    return None

# Options to consider
options = ["PK", "The Great Gatsby", "Philadelphia", "Casper Meets Wendy"]

# Find the most similar movie
similar_movie = find_similar_movie(options, target_genres)
print(similar_movie)