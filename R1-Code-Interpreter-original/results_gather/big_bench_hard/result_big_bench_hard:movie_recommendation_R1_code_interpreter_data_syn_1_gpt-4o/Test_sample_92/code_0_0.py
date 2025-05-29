# Define the movies and their genres/themes
movies = {
    "Death to Smoochy": ["dark comedy", "satire"],
    "The Artist": ["silent film", "romance", "drama"],
    "LA Confidential": ["neo-noir", "crime", "drama"],
    "New York": ["drama", "historical"],
    "I Love You": ["romance", "drama"]
}

# Define the target genres/themes
target_genres = ["drama", "thriller", "sci-fi", "fantasy", "psychological", "satire"]

# Find the movie with the most matching genres/themes
best_match = None
max_matches = 0

for movie, genres in movies.items():
    matches = len(set(genres) & set(target_genres))
    if matches > max_matches:
        best_match = movie
        max_matches = matches

print(best_match)