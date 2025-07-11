def find_hook_flash_champions():
    """
    This function determines which champions from a given list can perform a
    hook-flash to extend their ability's range.

    The specific mechanic is defined as: initiating the hook ability in one direction
    and then using the Flash summoner spell in the opposite direction during the
    cast animation, causing the hook to travel a longer distance.

    Based on game mechanics in Season 14 (2024):
    - Blitzcrank's Q is fired instantly; Flash doesn't extend it.
    - Nautilus's Q is fired instantly; Flash doesn't extend it.
    - Pyke's Q can be used with Flash, but it's for repositioning, not extending the hook's projectile length in this manner.
    - Thresh can Flash during his Q's wind-up animation, which moves the hook's
      origin point and effectively extends its range. This is the correct interaction.
    """
    capable_champions = ["Thresh"]

    # Format the result as a comma-separated string
    result = ",".join(capable_champions)

    print(result)

find_hook_flash_champions()