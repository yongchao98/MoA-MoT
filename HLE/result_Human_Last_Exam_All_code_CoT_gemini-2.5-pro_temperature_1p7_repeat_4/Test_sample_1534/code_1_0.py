def get_adjacent_brain_regions():
    """
    This function provides the names of brain regions adjacent posteriorly
    to the palliovisceral lobe in a dwarf cuttlefish.
    This information is based on established neuroanatomical knowledge.
    """
    
    # The palliovisceral lobe is a large, unpaired lobe in the posterior
    # subesophageal mass of the cuttlefish brain. It controls the viscera,
    # chromatophores, and the funnel/mantle for jet propulsion.
    
    # The lobes located directly posterior to it are primarily involved
    # in controlling the fins for locomotion and stabilization.
    
    posterior_adjacent_lobes = [
        "Posterior pedal lobe (often called the fin lobes)",
        "Magnocellular lobe"
    ]
    
    # Print the result in a formatted way.
    print("The brain regions that lie adjacent posteriorly to the palliovisceral lobe in the dwarf cuttlefish are:")
    for region in posterior_adjacent_lobes:
        print(f"- {region}")

# Execute the function to get the answer.
get_adjacent_brain_regions()