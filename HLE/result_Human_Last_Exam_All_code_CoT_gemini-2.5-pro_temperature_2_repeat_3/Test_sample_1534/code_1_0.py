import json

def find_adjacent_brain_regions():
    """
    This function simulates a query to a neuroanatomical database
    to find brain regions adjacent to a specified lobe in a specific direction.
    """

    # A simplified model of the dwarf cuttlefish brain's spatial organization.
    # The relationships are based on established cephalopod neuroanatomy.
    cuttlefish_brain_map = {
        "palliovisceral_lobe": {
            "description": "Part of the posterior subesophageal mass, controlling visceral organs, gills, and mantle.",
            "adjacent": {
                "anterior": ["brachial_lobe", "pedal_lobe"],
                "posterior": ["inferior_buccal_lobe", "posterior_basal_lobe"],
                "dorsal": ["optic_lobes"],
                "ventral": ["static_lobes"]
            }
        },
        "inferior_buccal_lobe": {
            "description": "Involved in feeding and control of the buccal mass."
        },
        "posterior_basal_lobe": {
            "description": "Involved in motor control, particularly of the head and arms."
        }
    }

    target_lobe = "palliovisceral_lobe"
    direction = "posterior"

    try:
        # Retrieve the list of lobes from the 'database'
        adjacent_lobes = cuttlefish_brain_map[target_lobe]["adjacent"][direction]

        # Format names for clear output
        target_lobe_formatted = target_lobe.replace('_', ' ').title()
        adjacent_lobes_formatted = [name.replace('_', ' ').title() for name in adjacent_lobes]

        print(f"The brain regions that lie adjacent posteriorly to the {target_lobe_formatted} are:")
        for lobe in adjacent_lobes_formatted:
            print(f"- {lobe}")

    except KeyError:
        print(f"Could not find information for the '{target_lobe}' in the '{direction}' direction.")

if __name__ == '__main__':
    find_adjacent_brain_regions()