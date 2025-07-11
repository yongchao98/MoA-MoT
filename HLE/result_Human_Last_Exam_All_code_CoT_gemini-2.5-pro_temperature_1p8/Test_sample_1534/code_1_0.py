def get_adjacent_brain_regions():
    """
    This function retrieves information about brain regions adjacent to the
    palliovisceral lobe in a dwarf cuttlefish based on known cephalopod neuroanatomy.
    """
    # The palliovisceral lobe is a major component of the posterior suboesophageal mass.
    # The regions lying directly posterior to it control the fins and posterior chromatophores.
    posterior_regions = {
        "Fin Lobes": "Control the movement of the fins.",
        "Posterior Chromatophore Lobes": "Control the color-changing skin cells (chromatophores) on the posterior part of the body."
    }
    return posterior_regions

def main():
    """
    Main function to execute the script.
    """
    print("Based on cephalopod neuroanatomy, the brain regions that lie adjacent posteriorly to the palliovisceral lobe in the dwarf cuttlefish are:")
    
    regions = get_adjacent_brain_regions()
    
    for region, description in regions.items():
        print(f"\n- Region: {region}")
        print(f"  Function: {description}")

if __name__ == "__main__":
    main()