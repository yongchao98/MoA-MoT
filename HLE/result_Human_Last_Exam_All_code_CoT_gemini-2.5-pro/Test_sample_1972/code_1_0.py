def find_borges_reference():
    """
    This script identifies and prints the novel and author
    Jorge Luis Borges was referring to with the given descriptions.
    """
    author = "Juan Carlos Onetti"
    novel = "El astillero (The Shipyard)"
    quote = "the intensity of a tiger and the variety that a chess duel can achieve"
    description = "a continuator (and simplifier) of Faulkner"

    print(f"The author Jorge Luis Borges described as '{description}' is {author}.")
    print(f"The novel he praised for having \"{quote}\" is '{novel}'.")

if __name__ == "__main__":
    find_borges_reference()