def count_triple_crown_winners():
    """
    This function lists the MLB Triple Crown winners,
    constructs an equation to count them, and prints the total.
    """
    # A list of unique players who have won the MLB Triple Crown.
    # Note: Rogers Hornsby and Ted Williams each won twice, but are counted as one winner each.
    winners = [
        "Paul Hines",
        "Tip O'Neill",
        "Hugh Duffy",
        "Nap Lajoie",
        "Ty Cobb",
        "Heinie Zimmerman",
        "Rogers Hornsby",
        "Chuck Klein",
        "Jimmie Foxx",
        "Lou Gehrig",
        "Joe Medwick",
        "Ted Williams",
        "Mickey Mantle",
        "Frank Robinson",
        "Carl Yastrzemski",
        "Miguel Cabrera"
    ]

    # The total number of winners is the length of the list.
    total_winners = len(winners)

    # Create a string representing the equation 1 + 1 + ... = total
    equation_string = " + ".join(['1'] * total_winners)

    print("The equation for counting the number of MLB Triple Crown winners is based on adding 1 for each unique player:")
    print(f"{equation_string} = {total_winners}")
    print(f"\nThere are {total_winners} MLB Triple Crown winners.")

# Run the function
count_triple_crown_winners()