def solve_heraldry_puzzle():
    """
    Identifies the city with a coat of arms similar to the one provided.
    """
    city_name = "Ananyiv, Ukraine"
    
    # Key elements from the image's coat of arms
    shield_division = "Horizontally divided (per fess)"
    top_field_color = "Blue (azure)"
    top_field_charge = "Gold crescent moon"
    bottom_field_color = "Gold (or)"
    bottom_field_charge = "Black crab"
    
    print(f"The city with a coat of arms that is very similar to the one in the image is: {city_name}")
    print("\nHere is a breakdown of the comparison:")
    print(f"1. Shield Division: The shield in the image is {shield_division}, just like the coat of arms of {city_name}.")
    print(f"2. Top Half: It features a {top_field_charge} on a {top_field_color} background, matching {city_name}'s coat of arms.")
    print(f"3. Bottom Half: It features a {bottom_field_charge} on a {bottom_field_color} background, which is also the main feature of {city_name}'s coat of arms.")
    print("\nThe crown and wreath are artistic embellishments, but the core heraldic design on the shield is a clear match.")

solve_heraldry_puzzle()