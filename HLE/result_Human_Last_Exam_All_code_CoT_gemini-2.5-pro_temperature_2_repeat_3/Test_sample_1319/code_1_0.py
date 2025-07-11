def find_tallest_iliad_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    It stores the mountains and their elevations, finds the tallest among them (excluding Olympus),
    and prints the result along with a comparison.
    """
    # A dictionary of mountains mentioned in the Iliad and their elevations in meters.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ossa": 1978,
        "Mount Ida (Turkey)": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611,
    }

    # Remove Mount Olympus for the comparison
    relevant_mountains = mountains.copy()
    if "Mount Olympus" in relevant_mountains:
        del relevant_mountains["Mount Olympus"]

    # Find the tallest mountain in the remaining list
    tallest_mountain_name = max(relevant_mountains, key=relevant_mountains.get)
    
    # Sort the mountains by height for the final equation
    sorted_mountains = sorted(relevant_mountains.items(), key=lambda item: item[1], reverse=True)

    print(f"After Mount Olympus, the tallest historical mountain mentioned in the Iliad is: {tallest_mountain_name}\n")
    print("Here is the equation comparing the mountain heights:")

    # Build and print the comparison string, showing each number
    equation_parts = []
    for name, height in sorted_mountains:
        equation_parts.append(f"{name} ({height}m)")
    
    final_equation = " > ".join(equation_parts)
    print(final_equation)

if __name__ == '__main__':
    find_tallest_iliad_mountain()