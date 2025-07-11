def count_milton_saints():
    """
    Calculates and prints the number of historical saints mentioned by name
    in John Milton's "Paradise Lost".

    The list of saints is pre-compiled from a textual analysis of the poem,
    focusing on biblical figures (patriarchs, prophets, angels) who are
    venerated as saints in major Christian traditions.
    """
    
    # Based on a literary analysis of "Paradise Lost", these are the
    # historical/biblical figures venerated as saints mentioned by name.
    named_saints = [
        "Adam",      # First man, venerated as a saint
        "Eve",       # First woman, venerated as a saint
        "Abel",      # Son of Adam and Eve, first martyr
        "Enoch",     # Patriarch who "walked with God"
        "Noah",      # Patriarch who built the Ark
        "Abraham",   # Patriarch of the Abrahamic religions
        "Moses",     # Prophet who led the Israelites
        "Joshua",    # Successor to Moses
        "David",     # King of Israel and psalmist
        "Michael",   # Archangel and commander of God's army
        "Gabriel",   # Archangel and messenger
        "Raphael",   # Archangel who guides Adam
        "Uriel"      # Archangel, Regent of the Sun
    ]

    # Calculate the total number of saints found
    count = len(named_saints)

    # Format the output to show the individual names contributing to the sum
    equation_str = " + ".join(["1" for _ in named_saints])
    names_str = ", ".join(named_saints)
    
    print(f"The number of historical saints mentioned by name in 'Paradise Lost' is {count}.")
    print(f"They are: {names_str}.")
    print("\nThis can be represented as the following equation:")
    print(f"{equation_str} = {count}")


if __name__ == "__main__":
    count_milton_saints()
<<<13>>>