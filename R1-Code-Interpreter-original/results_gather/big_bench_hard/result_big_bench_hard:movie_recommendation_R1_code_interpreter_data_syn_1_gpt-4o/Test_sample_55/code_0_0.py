# Define a dictionary with movie options and their characteristics
movies = {
    "Hawks and Sparrows": {"genre": "comedy", "popularity": "low"},
    "The Scent of Green Papaya": {"genre": "drama", "popularity": "medium"},
    "Schindler's List": {"genre": "historical drama", "popularity": "high"},
    "Le Maître d'école": {"genre": "comedy", "popularity": "low"}
}

# Define the characteristics of the given movies
given_movies = {
    "Terminator 2": {"genre": "action", "popularity": "high"},
    "Pulp Fiction": {"genre": "crime", "popularity": "high"},
    "Aladdin": {"genre": "animated musical", "popularity": "high"},
    "The Lion King": {"genre": "animated musical", "popularity": "high"}
}

# Find the most similar movie based on popularity
most_similar = max(movies, key=lambda x: movies[x]["popularity"] == "high")

print(most_similar)