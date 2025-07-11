def find_runestone_id():
    """
    Identifies the Ingvar runestone from an inscription fragment.
    """
    # A small database of Ingvar runestone inscriptions (fragments)
    ingvar_runestones = {
        "Sö 9": "þuriʀ : risti : runaʀ : þasi : eftiʀ : tuka : sun : sin : auk : hlaþkair",
        "Sö 179": "tula : lit : raisa : stain : þinsa| |at : sun : sin : haralt : bruþur : inkuars",
        "U 654": "alrikr : raisti : stain : sun : siriþaʀ : at : sin : faþur : spak : ... miþ : ikuari",
        "U 778": "anduitr : auk : kiti : auk : kari : auk : blisi : þaiʀ : raistu : stain : þina : at : kunar : faþur sin : ... austr : miþ : ikuari"
    }

    # The inscription fragment identified from the image
    image_fragment = "þuriʀ : risti : runaʀ"

    found_id = None
    for stone_id, inscription in ingvar_runestones.items():
        if image_fragment in inscription:
            found_id = stone_id
            break
    
    if found_id:
        print(f"The inscription fragment belongs to the Ingvar runestone with the ID: {found_id}")
    else:
        print("Could not identify the runestone from the provided fragment.")

find_runestone_id()