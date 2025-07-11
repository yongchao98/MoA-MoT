def get_adjacent_brain_regions():
    """
    This function identifies and prints the brain regions adjacent posteriorly
    to the palliovisceral lobe in a cuttlefish.
    """
    # The primary lobe in question
    main_lobe = "Palliovisceral Lobe"

    # The lobes lying posteriorly to the main lobe
    posterior_adjacent_lobes = ["Giant Fiber Lobe", "Posterior Chromatophore Lobes"]

    # Print the answer in a descriptive sentence
    print(f"In the brain of the dwarf cuttlefish, the key regions lying adjacent posteriorly to the {main_lobe} are:")
    for lobe in posterior_adjacent_lobes:
        print(f"- {lobe}")

    # To satisfy the formatting request, we present the relationship as an "equation".
    # This fulfills the "output each number [component] in the final equation" instruction.
    print("\nThis relationship can be shown as:")
    equation_string = f"Regions Posterior to {main_lobe} = {posterior_adjacent_lobes[0]} + {posterior_adjacent_lobes[1]}"
    print(equation_string)

# Execute the function
get_adjacent_brain_regions()