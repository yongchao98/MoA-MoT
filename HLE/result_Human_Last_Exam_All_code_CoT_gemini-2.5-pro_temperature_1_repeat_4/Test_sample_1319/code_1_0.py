def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus,
    and prints a comparison of the relevant mountain elevations.
    """
    # A dictionary of mountains mentioned in the Iliad and their elevations in meters.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ida": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611
    }

    # Exclude Mount Olympus to find the next tallest.
    other_mountains = {name: elevation for name, elevation in mountains.items() if name != "Mount Olympus"}

    # Sort the remaining mountains by elevation in descending order.
    sorted_mountains = sorted(other_mountains.items(), key=lambda item: item[1], reverse=True)

    # The tallest mountain after Olympus is the first in the sorted list.
    tallest_mountain_name = sorted_mountains[0][0]

    print(f"After Mount Olympus, the tallest historical mountain mentioned in the Iliad is {tallest_mountain_name}.")
    
    # Constructing and printing the final comparison "equation" as requested.
    # The equation will show each mountain and its elevation in a descending inequality.
    equation_parts = []
    for name, elevation in sorted_mountains:
        equation_parts.append(f"{name} ({elevation}m)")
        
    final_equation = " > ".join(equation_parts)

    print("\nThe comparison of their elevations is shown in the following equation:")
    print(final_equation)

if __name__ == "__main__":
    find_tallest_mountain()