def find_hook_flash_champions():
    """
    Identifies which champions from a given list can perform a range-extending hook-flash.

    The specific mechanic is:
    1. Cast the hook ability in one direction.
    2. During the cast animation, use Flash in the opposite direction.
    3. The hook originates from the new location but travels in the original direction, extending its effective range.
    """
    champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    
    # Based on game mechanics as of Season 14 (2024):
    # - Blitzcrank's Q is a projectile; flashing after cast does nothing.
    # - Nautilus's Q is a projectile; flashing after cast does nothing.
    # - Pyke's Q is a projectile; flashing after cast does nothing.
    # - Thresh's Q direction is locked on press, but its origin is determined at the end of the cast time,
    #   allowing Flash to change the origin and extend the hook's range.
    
    capable_champions = ["Thresh"]
    
    if not capable_champions:
        result = "None"
    else:
        result = ", ".join(capable_champions)
        
    print(result)

find_hook_flash_champions()