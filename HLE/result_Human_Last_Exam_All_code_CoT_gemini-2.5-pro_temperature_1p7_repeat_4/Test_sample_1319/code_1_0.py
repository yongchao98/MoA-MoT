import sys

def find_tallest_mountain():
    """
    This script identifies the tallest historical mountain mentioned in the Iliad,
    excluding Mount Olympus, by comparing the elevations of known candidates.
    """
    # A dictionary of historical mountains from the Iliad (excluding Olympus) 
    # and their approximate highest point in meters.
    mountains = {
        "Mount Ida": 1774,
        "Samothrace (Mount Fengari)": 1611,
        "Mount Athos": 2033,
        "Mount Ossa": 1978,
        "Mount Pelion": 1624
    }

    # Find the name and height of the tallest mountain in the dictionary
    tallest_mountain_name = ""
    max_height = -1

    for name, height in mountains.items():
        if height > max_height:
            max_height = height
            tallest_mountain_name = name

    # Print the list of mountains and their heights for context
    print("The historical mountains mentioned in the Iliad (besides Olympus) and their approximate heights are:")
    for name, height in mountains.items():
        print(f"- {name}: {height} meters")

    # Get the heights for the final "equation" output
    # The map(str, ...) converts each integer height to a string
    # The ", ".join(...) joins them together with a comma and space
    heights_str = ", ".join(map(str, mountains.values()))
    
    print("\nTo determine the tallest, we find the maximum of these heights:")
    # This print statement fulfills the requirement to output each number in the final equation.
    print(f"max({heights_str}) = {max_height}")
    
    print(f"\nTherefore, the tallest historical mountain is {tallest_mountain_name}.")

find_tallest_mountain()