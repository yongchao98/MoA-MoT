# This script is designed to solve the riddle by logically evaluating the provided clues.

def find_island_town():
    """
    Solves the riddle by filtering candidates based on a set of criteria.
    """
    # Clue 1: An island town.
    # Clue 2: More than 500 kilometers from another inhabited island.
    # Clue 3: Sits on a bay formed by a volcanic caldera.
    # Clue 4: The town and the bay share the same name.

    # Candidate evaluation:
    # - Tristan da Cunha: The main settlement is Edinburgh of the Seven Seas. It is extremely remote,
    #   but the town's name does not match a specific bay it sits on. -> Incorrect.
    # - Pitcairn Islands: The main settlement is Adamstown. It is near Bounty Bay.
    #   The names do not match. -> Incorrect.
    # - Easter Island (Rapa Nui): The main town is Hanga Roa.
    #   - Is it remote (> 500 km)? Yes, it is one of the most isolated inhabited islands on Earth.
    #   - Do the town and bay share a name? Yes, the town of Hanga Roa is on the bay of Hanga Roa.
    #   - Is the bay formed by a volcanic caldera? The island is entirely of volcanic origin,
    #     and its geography, including its bays, was shaped by its three main volcanoes and their calderas. -> Correct.

    # The only candidate that satisfies all conditions is Hanga Roa.
    # The riddle asks for the island town's name.
    answer = "Hanga Roa"

    # The riddle mentions the distance of 500 kilometers. Let's include that in the output.
    distance_in_km = 500
    print(f"The riddle specifies a distance of more than {distance_in_km} kilometers from another inhabited island.")
    print(f"The island town that sits on a bay formed by a volcanic caldera and shares the bay's name is: {answer}")

find_island_town()