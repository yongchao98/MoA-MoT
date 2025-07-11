import sys

def find_movie_from_trivia(clue_key):
    """
    Looks up a movie title from a predefined trivia knowledge base.
    """
    # A small knowledge base mapping trivia clues to movie titles.
    trivia_knowledge_base = {
        "bus_fate_thank_you": {
            "question": "In which Oscar-nominated film does a character seal his fate by saying \"Thank you\" while boarding a bus?",
            "film": "Midnight Cowboy",
            "character": "Ratso Rizzo",
            "context": "At the end of the film, a very ill Ratso Rizzo boards a bus with Joe Buck to go to Florida. This final journey is where he passes away, effectively sealing his fate."
        }
    }

    # Retrieve the trivia item from the knowledge base.
    trivia_item = trivia_knowledge_base.get(clue_key)

    if trivia_item:
        print(f"Query: {trivia_item['question']}")
        print("-" * 20)
        print(f"Answer: The film is '{trivia_item['film']}'.")
        print(f"Context: In the film, the character {trivia_item['character']} dies on a bus ride to Florida, which he takes in a last-ditch effort to find a better life.")
    else:
        print(f"Sorry, the clue '{clue_key}' was not found in the knowledge base.", file=sys.stderr)

if __name__ == "__main__":
    # The key corresponding to the user's question.
    key_for_clue = "bus_fate_thank_you"
    find_movie_from_trivia(key_for_clue)
