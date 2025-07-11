def find_posterior_brain_regions():
    """
    This function identifies and prints the brain regions located posteriorly
    to the palliovisceral lobe in a cuttlefish.
    
    The palliovisceral lobe is a part of the posterior subesophageal mass
    of the cephalopod brain. The structures lying adjacent and posterior to it
    are critical for locomotion and escape responses.
    """
    
    # Define the brain regions posterior to the palliovisceral lobe
    posterior_regions = ["Fin lobes", "Giant fiber lobe"]
    
    print("The brain regions that lie adjacent posteriorly to the palliovisceral lobe in the dwarf cuttlefish are:")
    
    # Print each region
    for region in posterior_regions:
        print(f"- {region}")

# Execute the function to display the answer
find_posterior_brain_regions()