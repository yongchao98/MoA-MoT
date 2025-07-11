import collections

def count_phosphorus_colors():
    """
    This function identifies and counts the number of distinct colors
    observed in the common pure allotropes of phosphorus.
    """
    # The primary, distinctively colored allotropes of phosphorus are:
    # 1. White Phosphorus (can also be yellowish)
    # 2. Red Phosphorus
    # 3. Violet Phosphorus (also known as Hittorf's metallic phosphorus)
    # 4. Black Phosphorus (the most thermodynamically stable form)
    
    allotropes = {
        "White/Yellow": 1,
        "Red": 1,
        "Violet": 1,
        "Black": 1
    }

    print("The common allotropes of phosphorus and their colors are:")
    for color in allotropes:
        print(f"- {color} phosphorus")
    
    print("\nTo count the total number of colors, we can form an equation:")
    
    # Create the equation string from the numbers
    numbers = list(allotropes.values())
    total = sum(numbers)
    
    equation_string = " + ".join(map(str, numbers))
    
    print(f"{equation_string} = {total}")

count_phosphorus_colors()