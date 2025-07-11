def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, after Mount Olympus.
    """
    # Step 1 & 2: Define the mountains mentioned in the Iliad and their elevations in meters.
    mountains = {
        "Mount Olympus": 2918,
        "Mount Athos": 2033,
        "Mount Ida (Turkey)": 1774,
        "Mount Samothrace": 1611,
    }

    print("Mountains mentioned in the Iliad (with elevations):")
    for name, height in mountains.items():
        print(f"- {name}: {height}m")
    print("-" * 20)

    # Step 3: Exclude Mount Olympus from the comparison.
    contenders = {name: height for name, height in mountains.items() if name != "Mount Olympus"}

    # Step 4: Find the tallest among the remaining mountains.
    if not contenders:
        print("No other mountains to compare.")
        return

    # Sort contenders by height in descending order to create the "equation"
    sorted_contenders = sorted(contenders.items(), key=lambda item: item[1], reverse=True)
    
    tallest_mountain_name = sorted_contenders[0][0]
    tallest_mountain_height = sorted_contenders[0][1]

    print("After excluding Mount Olympus, the comparison is:")
    
    # Building the final equation string as requested
    equation_parts = []
    for name, height in sorted_contenders:
        equation_parts.append(f"{name} ({height}m)")
    
    final_equation = " > ".join(equation_parts)
    print(final_equation)
    print("-" * 20)
    
    print(f"The tallest historical mountain mentioned in the Iliad after Mount Olympus is {tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")

find_tallest_mountain()