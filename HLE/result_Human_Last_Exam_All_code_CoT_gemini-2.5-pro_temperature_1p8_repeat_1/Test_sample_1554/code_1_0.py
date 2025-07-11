def find_philosopher():
    """
    This function identifies and prints the name of the Lviv-Warsaw School philosopher
    who distinguished between symptoms and signals.
    """
    philosopher = "Kazimierz Ajdukiewicz"
    concept = "hallmark (or sign)"
    decomposition = {
        "symptoms": "natural signs with a causal link to their object",
        "signals": "conventional signs based on rules or agreement"
    }

    print(f"The philosopher from the Lviv-Warsaw School who decomposes the notion of a '{concept}' is {philosopher}.")
    print(f"He broke it down into:")
    print(f"- Symptoms: {decomposition['symptoms']}")
    print(f"- Signals: {decomposition['signals']}")

# Execute the function to display the answer
find_philosopher()