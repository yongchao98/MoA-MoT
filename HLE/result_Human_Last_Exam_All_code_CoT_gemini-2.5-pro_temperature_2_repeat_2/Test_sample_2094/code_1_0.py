def find_hook_flash_champions():
    """
    This script determines which champions from a predefined list can perform a hook-flash.

    The "hook-flash" is a well-known League of Legends mechanic. It involves casting a hook
    ability and then using the Flash summoner spell during the ability's startup animation.
    This action changes the hook's origin point to the post-Flash location, which can
    extend the ability's range or change its angle to catch opponents off-guard.

    All four champions listed—Blitzcrank, Nautilus, Pyke, and Thresh—are famous for
    their ability to perform this technique with their Q abilities.
    """
    
    # List of champions to check
    champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]

    # In League of Legends (as of Season 14), all the listed champions can use
    # Flash during the cast animation of their hook ability (Q) to extend its range.
    capable_champions = champions

    # Prepare the final answer string
    final_answer = ",".join(capable_champions)

    print("The following champions can perform a hook-flash to extend their hook's range:")
    print(final_answer)

find_hook_flash_champions()