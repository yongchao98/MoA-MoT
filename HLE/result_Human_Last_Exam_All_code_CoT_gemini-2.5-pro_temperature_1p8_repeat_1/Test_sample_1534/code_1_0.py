def get_posterior_regions():
    """
    This function identifies and prints the brain structures located
    posteriorly to the palliovisceral lobe in a cuttlefish.
    """
    # The palliovisceral lobe is the most posterior major lobe of the central brain mass.
    # Therefore, the structures lying adjacent and posterior to it are primarily the
    # nerve-originating regions within it and the major nerve trunks leaving it.

    # The sub-region within the posterior palliovisceral lobe that gives rise to the giant axons.
    region1 = "Giant Fiber Lobes"

    # The major nerve trunks exiting the brain posteriorly.
    region2 = "Posterior Pallial Nerves"
    region3 = "Visceral Nerves"

    print("In the dwarf cuttlefish brain, the structures lying adjacent posteriorly to the main body of the palliovisceral lobe are:")
    print(f"- The {region1}, which are located in the posterior part of the lobe complex and give rise to the giant axons.")
    print(f"- The large nerve trunks that exit the brain, primarily the {region2} and the {region3}.")

if __name__ == '__main__':
    get_posterior_regions()