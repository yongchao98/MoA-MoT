def solve_coat_of_arms_riddle():
    """
    This script analyzes the provided coat of arms and identifies the city
    with a very similar emblem.
    """
    # Step 1: Analyze the key features of the coat of arms from the image.
    shield_division = "Horizontally divided (per fess)"
    top_half_color = "Blue (azure)"
    top_half_charge = "Crescent moon (gold)"
    bottom_half_color = "Gold (or)"
    bottom_half_charge = "Crab (black)"

    print("Analysis of the Coat of Arms:")
    print(f"- Shield Division: {shield_division}")
    print(f"- Top Half: A {top_half_charge} on a {top_half_color} field.")
    print(f"- Bottom Half: A {bottom_half_charge} on a {bottom_half_color} field.")
    print("-" * 20)

    # Step 2: Identify the corresponding real-world emblem.
    print("Identification:")
    print("The combination of these elements is very specific.")
    print("It matches the blazon (formal description) of the coat of arms for Ilfov County in Romania.")
    print("Historically, this same design was the coat of arms for Constanța County, Romania.")
    print("-" * 20)

    # Step 3: Connect the emblem to a specific city as requested.
    print("Connecting to a City:")
    print("The prompt asks for a city, not a county.")
    print("The symbol of a crab is highly relevant to a coastal or port city.")
    print("Constanța is a major port city on the Black Sea and was the historical user of this county-level emblem.")
    print("The symbolism is much more strongly tied to Constanța than to the landlocked Ilfov County.")
    print("-" * 20)

    # Step 4: Conclude and print the answer.
    final_city = "Constanța"
    print(f"Conclusion: The city with a coat of arms very similar to the one in the image is Constanța.")
    print("\nFinal Answer:")
    print(final_city)

# Run the solver
solve_coat_of_arms_riddle()