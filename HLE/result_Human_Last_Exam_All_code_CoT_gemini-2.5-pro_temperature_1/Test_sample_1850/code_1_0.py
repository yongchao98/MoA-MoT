import collections

def count_milton_saints():
    """
    Identifies and counts the historical saints mentioned by name in Milton's "Paradise Lost".
    """
    # Based on a literary analysis of "Paradise Lost", these are the historical
    # figures treated as saints who are mentioned by their proper name.
    # Most are found in the vision shown to Adam in Books XI and XII.
    saints = collections.OrderedDict([
        ("Adam", "Central figure and first man"),
        ("Eve", "Central figure and first woman"),
        ("Abel", "Righteous son of Adam, slain by Cain (Book XI)"),
        ("Enoch", "Patriarch who was taken by God (Book XI)"),
        ("Noah", "Builder of the Ark (Book XI)"),
        ("Abraham", "Patriarch of the Hebrew race (Book XII)"),
        ("Isaac", "Son of Abraham (Book XII)"),
        ("Jacob", "Son of Isaac (Book XII)"),
        ("Joseph", "Son of Jacob (Book XII)"),
        ("Moses", "Leader of the Exodus (Book XII)"),
        ("Joshua", "Successor to Moses (Book XII)"),
        ("David", "King of Israel (Book XII)"),
        ("Solomon", "Son of David, builder of the Temple (Book XII)")
    ])

    count = len(saints)

    print("The historical saints mentioned by name in Milton's 'Paradise Lost' are:")
    for name, description in saints.items():
        print(f"- {name}")

    # Create the equation string as requested
    equation_parts = ["1"] * count
    equation_str = " + ".join(equation_parts)

    print(f"\nThis gives a total of {count} named historical saints.")
    print("\nThe final calculation is:")
    print(f"{equation_str} = {count}")

if __name__ == "__main__":
    count_milton_saints()