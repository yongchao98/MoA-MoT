def find_unmentioned_odyssey_elements():
    """
    This script determines which chemical elements, known to the ancient Greeks
    at the time of the Odyssey's composition, are not mentioned in the epic.
    """
    # Step 1: Define the set of elements known in antiquity (c. 8th century BCE),
    # excluding zinc and antimony as per the problem description.
    known_elements_in_antiquity = {"Carbon", "Sulfur", "Iron", "Copper", "Silver", "Tin", "Gold", "Mercury", "Lead"}

    # Step 2: Define the set of elements from the above list that are mentioned
    # in Homer's Odyssey.
    # 'khalkos' (χαλκός) refers to copper/bronze, 'kassiteros' (κασσίτερος) to tin,
    # 'sideros' (σίδηρος) to iron, 'theion' (θεῖον) to sulfur, and 'anthrax' (ἄνθραξ) to charcoal/carbon.
    # Gold and silver are mentioned extensively.
    elements_mentioned_in_odyssey = {"Carbon", "Sulfur", "Iron", "Copper", "Silver", "Tin", "Gold"}

    # Step 3: Find the difference between the two sets to isolate the unmentioned elements.
    unmentioned_elements = known_elements_in_antiquity.difference(elements_mentioned_in_odyssey)

    # Convert to a sorted list for consistent and clean output.
    final_list = sorted(list(unmentioned_elements))

    print("Of the elements known in their elemental form during the time the Odyssey was composed, the following are not mentioned in the poem:")
    for element in final_list:
        print(f"- {element}")

if __name__ == "__main__":
    find_unmentioned_odyssey_elements()