def find_mutualists():
    """
    This function identifies the indices of mutualists of Asclepias fascicularis
    from a predefined list and prints them.

    A mutualistic relationship is one where both organisms benefit.
    For Asclepias fascicularis (milkweed), this primarily involves pollination.

    - Adults (Pollinators/Defenders):
      1) Danaus plexipus (Monarch butterfly) - Pollinator
      2) Megachile frigidus (Leafcutter bee) - Pollinator
      3) Formica rufa (Ant) - Nectar feeder, can act as a defender
      4) Sphex ichneumoneus (Digger wasp) - Pollinator
      5) Pepsis thisbe (Tarantula hawk wasp) - Pollinator
      6) Megachile ericetorum (Leafcutter bee) - Pollinator
    These are all mutualists.

    - Larvae:
      7) Danaus plexipus (Monarch larva) - Herbivore (eats the plant)
      8-12) Other larvae - Do not directly interact with the plant in a beneficial way.
    None of the larvae are mutualists.
    """
    mutualist_indices = [1, 2, 3, 4, 5, 6]
    
    # Convert numbers to strings to join them
    result_string = ",".join(map(str, mutualist_indices))
    
    print(result_string)

find_mutualists()