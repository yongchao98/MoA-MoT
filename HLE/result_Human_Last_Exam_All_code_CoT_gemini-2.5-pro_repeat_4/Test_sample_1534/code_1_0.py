def get_cuttlefish_anatomy():
    """
    This function provides information on the brain regions adjacent posteriorly
    to the palliovisceral lobe in the dwarf cuttlefish.

    The palliovisceral lobe is one of the most posterior lobes in the main
    brain mass of a cuttlefish. Therefore, the structures lying immediately
    posterior to it are not other lobes, but the major nerve trunks that
    emerge from it to innervate the rest of the body.
    """
    region1 = "pallial nerves"
    region2 = "visceral nerves"

    explanation = (
        f"The primary structures lying adjacent posteriorly to the palliovisceral lobe "
        f"in the dwarf cuttlefish are the large nerve trunks that originate from it. "
        f"These are the {region1} (which contain the giant axons controlling the mantle for jet propulsion) "
        f"and the {region2} (which control the internal organs)."
    )

    print(explanation)

get_cuttlefish_anatomy()