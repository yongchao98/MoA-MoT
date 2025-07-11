def find_adjacent_brain_regions():
    """
    This function stores and prints the known brain regions adjacent
    posteriorly to the palliovisceral lobe in the dwarf cuttlefish.
    """
    # Anatomical data based on cephalopod neuroanatomy
    target_lobe = "palliovisceral lobe"
    adjacent_posterior_regions = [
        "Posterior chromatophore lobes",
        "Giant fiber lobes"
    ]

    # Print the information in a clear statement
    print(f"The brain regions that lie adjacent posteriorly to the {target_lobe} in the dwarf cuttlefish are:")
    for region in adjacent_posterior_regions:
        print(f"- {region}")

# Execute the function to display the answer
find_adjacent_brain_regions()