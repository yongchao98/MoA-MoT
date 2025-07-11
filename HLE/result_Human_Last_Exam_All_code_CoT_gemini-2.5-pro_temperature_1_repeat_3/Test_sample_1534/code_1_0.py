def get_adjacent_brain_regions():
    """
    This function identifies and prints the brain regions located posteriorly
    to the palliovisceral lobe in a dwarf cuttlefish.
    """
    # The palliovisceral lobe is a large structure in the posterior subesophageal mass of the cephalopod brain.
    # It controls visceral organs and jet propulsion.
    # Anatomical studies show that motor control lobes for fins and chromatophores are situated just behind it.

    posteriorly_adjacent_lobes = {
        "Region 1": "Posterior chromatophore lobes",
        "Region 2": "Fin lobes"
    }

    print("In the dwarf cuttlefish, the brain regions lying adjacent posteriorly to the palliovisceral lobe are:")
    print(f"- {posteriorly_adjacent_lobes['Region 1']}")
    print(f"- {posteriorly_adjacent_lobes['Region 2']}")
    print("\nThese lobes are part of the posterior subesophageal mass and are primarily involved in motor control for coloration and fin movement.")

# Execute the function to display the information.
get_adjacent_brain_regions()