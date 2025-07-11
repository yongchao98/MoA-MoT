def find_adjacent_brain_regions():
    """
    This function identifies and prints the brain regions adjacent posteriorly
    to the palliovisceral lobe in a cuttlefish.
    
    The brain anatomy of cephalopods is complex. The palliovisceral lobe is
    part of the posterior subesophageal mass, a region controlling the mantle,
    funnel, and internal organs.
    
    Based on scientific brain atlases for Sepia (cuttlefish), the brain regions
    located immediately posterior to the palliovisceral lobe are the
    posterior basal lobe and the median basal lobe.
    """
    
    target_lobe = "palliovisceral lobe"
    adjacent_posterior_regions = [
        "posterior basal lobe",
        "median basal lobe"
    ]

    print(f"The brain regions that lie adjacent posteriorly to the {target_lobe} are:")
    for region in adjacent_posterior_regions:
        print(f"- {region}")

# Execute the function to provide the answer.
find_adjacent_brain_regions()