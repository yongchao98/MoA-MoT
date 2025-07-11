def identify_coin_origin():
    """
    Identifies the origin of the provided French colonial coin.
    """
    origin = "New France (modern-day Quebec, Canada)"
    inscription = "COLONIES FRANCOISES"
    ruler = "Louis XV"
    mint = "Paris (Mint Mark 'A')"
    denomination = "9 Deniers"
    year = "1721-1722"

    print(f"The coin shown is a {denomination} from the reign of King {ruler}.")
    print(f"It was minted in {mint} during the years {year}.")
    print(f"The inscription '{inscription}' indicates it was created for circulation in the French colonies.")
    print("Specifically, this issue was intended for the North American territory known as:")
    print(f"Origin: {origin}")

identify_coin_origin()