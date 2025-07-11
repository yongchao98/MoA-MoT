def get_adjacent_brain_regions():
    """
    This function identifies and prints the brain regions located posteriorly
    to the palliovisceral lobe in a dwarf cuttlefish.
    The information is based on established cephalopod neuroanatomy.
    """
    palliovisceral_lobe_location = "part of the posterior subesophageal mass"
    posteriorly_adjacent_regions = [
        "First-order giant fiber lobes",
        "Second-order giant fiber lobes",
        "Third-order giant fiber lobes"
    ]

    print("The palliovisceral lobe is a major component of the subesophageal mass in the cuttlefish brain.")
    print("\nThe brain regions lying adjacent posteriorly to the palliovisceral lobe are:")
    for region in posteriorly_adjacent_regions:
        print(f"- {region}")

get_adjacent_brain_regions()