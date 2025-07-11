def find_mutualists():
    """
    This function identifies the mutualists of Asclepias fascicularis from a predefined list.

    A mutualist relationship benefits both organisms. In this case, the adult insects that
    pollinate the plant in exchange for nectar are considered mutualists.

    - Adults 1 (Danaus plexipus), 2 (Megachile frigidus), 4 (Sphex ichneumoneus),
      5 (Pepsis thisbe), and 6 (Megachile ericetorum) are all pollinators.
    - Adult 3 (Formica rufa) is an ant, generally considered a poor pollinator and sometimes
      protects plant pests like aphids. It is not a clear mutualist.
    - Larva 7 (Danaus plexipus) is a herbivore that eats the plant's leaves (antagonism).
    - Larvae 8, 9, 10, 11, and 12 do not directly interact with the living plant.

    Therefore, the mutualists are the adult pollinators.
    """
    mutualist_indices = [1, 2, 4, 5, 6]

    # Convert the list of numbers into a comma-separated string for output.
    output_string = ",".join(map(str, mutualist_indices))

    print(output_string)

find_mutualists()