def get_adjacent_brain_regions():
    """
    This function provides the names of brain regions adjacent posteriorly
    to the palliovisceral lobe in a cuttlefish.
    """
    # The palliovisceral lobe is a very posterior part of the subesophageal mass.
    # The regions located adjacent and posterior to it are primarily the fin lobes
    # and the posterior chromatophore lobes, which control the fins and skin patterning, respectively.
    adjacent_posterior_regions = ["posterior chromatophore lobes", "fin lobes"]

    print("The brain regions that lie adjacent posteriorly to the palliovisceral lobe in the dwarf cuttlefish are the:")
    for region in adjacent_posterior_regions:
        print(f"- {region}")

if __name__ == "__main__":
    get_adjacent_brain_regions()