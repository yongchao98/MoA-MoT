def find_painting_title():
    """
    Identifies a painting by matching its described features against a known database.
    """
    # Database of paintings with their key features (title, artist, century, keywords)
    paintings_db = {
        "The Mill": {
            "artist": "Rembrandt",
            "century": 17,
            "keywords": {"windmill", "dutch", "dark", "river"}
        },
        "The Windmill at Wijk bij Duurstede": {
            "artist": "Jacob van Ruisdael",
            "century": 17,
            "keywords": {"windmill", "dutch", "clouds", "river"}
        },
        "The Burning Mill": {
            "artist": "August Strindberg",
            "century": 19,
            "keywords": {"windmill", "fire", "storm", "dramatic"}
        },
        "Windmill at Brighton": {
            "artist": "John Constable",
            "century": 19,
            "keywords": {"windmill", "english", "landscape"}
        }
    }

    # Features identified from the image provided by the user
    image_features = {
        "century": 19,
        "keywords": {"windmill", "fire", "storm"}
    }

    best_match_title = "Unknown"
    highest_score = 0

    # Find the best match from the database
    for title, data in paintings_db.items():
        score = 0
        if data["century"] == image_features["century"]:
            # Give a higher weight to matching the century
            score += 2
        
        # Calculate keyword match score
        score += len(image_features["keywords"].intersection(data["keywords"]))

        if score > highest_score:
            highest_score = score
            best_match_title = title

    print(f"The title of the painting is: {best_match_title}")

find_painting_title()