import operator

def find_tallest_mountain():
    """
    Identifies the tallest historical mountain mentioned in the Iliad, besides Mount Olympus.
    """
    # Step 1, 2, 3: A dictionary of mountains from the Iliad (excluding Olympus)
    # and their corresponding real-world elevations in meters.
    # Mount Ida (Kaz Dağı) in Turkey, from where the gods watched the Trojan War.
    # Mount Athos in Greece, which Hera flew past.
    # Mount Samothrace (Fengari), where Poseidon sat to watch the war.
    # Mount Pelion in Greece, from where the wood for Achilles' spear came.
    mountains = {
        "Mount Athos": 2033,
        "Mount Ida": 1774,
        "Mount Pelion": 1624,
        "Mount Samothrace": 1611
    }

    # Step 4: Compare the elevations to find the tallest mountain.
    tallest_mountain_name = max(mountains.items(), key=operator.itemgetter(1))[0]
    tallest_mountain_height = mountains[tallest_mountain_name]

    print("Comparing the elevations of historical mountains from the Iliad (other than Olympus):")
    # "output each number in the final equation!"
    # The "equation" here is the comparison of the heights.
    for name, height in mountains.items():
        print(f"- {name}: {height} meters")

    print(f"\nThe tallest among these is {tallest_mountain_name} with an elevation of {tallest_mountain_height} meters.")


if __name__ == "__main__":
    find_tallest_mountain()