import sys

def identify_jazz_standard():
    """
    This function identifies the jazz standard based on its melodic and harmonic content.
    """
    title = "Half Nelson"
    composer = "Miles Davis"
    original_composition_title = "Lady Bird"
    original_composition_composer = "Tadd Dameron"

    print(f"Title: {title}")
    print(f"Composer: {composer}")
    print("\n# Note:")
    print("# This piece is a famous 'contrafact'.")
    print(f"# The melody is '{title}' by {composer}, but it is written over the chord changes of")
    print(f"# the jazz standard '{original_composition_title}' by {original_composition_composer}.")

if __name__ == "__main__":
    identify_jazz_standard()
