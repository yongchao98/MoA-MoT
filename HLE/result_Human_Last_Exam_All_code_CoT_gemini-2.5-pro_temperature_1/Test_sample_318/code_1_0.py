def find_synonym_indices():
    """
    This function identifies and prints the indices of species from the 1872 list
    that are considered synonyms in modern taxonomy.

    The synonym data is based on research from modern taxonomic databases. A name
    is considered a synonym if it's no longer the accepted scientific name for a species.
    This includes changes in genus, species epithet, or if the name was found to be
    a junior synonym of a previously described species.
    """

    # A list of the species names from the 1872 survey. The index in this list is (number - 1).
    species_list_1872 = [
        "Cimbex americana, var. Ulmi", # Valid species, variety is a synonym of the species.
        "Abia Kennicotti",             # Synonym of Abia americana
        "Acordulecera dorsalis",       # Valid
        "Ptenos texanus",              # Synonym of Sterictiphora texana
        "Ptenos niger",                # Synonym of Sterictiphora niger
        "Ptenos nigropectus",          # Synonym of Sterictiphora nigropectus
        "Hylotoma abdominalis",        # Synonym of Arge abdominalis
        "Hylotoma miniata",            # Synonym of Arge miniata
        "Hylotoma rubiginosa",         # Synonym of Arge rubiginosa
        "Nematus chloreus",            # Synonym of Pristiphora chlorea
        "Emphytus Bollii",             # Synonym of Ametastegia bollii
        "Hemichroa albidovariata",     # Valid
        "Hemichroa fraternalis",       # Synonym of Caulocampus fraterna
        "Selandria inaequidens",       # Synonym of Eriocampa inaequidens
        "Selandria albicollis",        # Synonym of Periclista albicollis
        "Macrophya excavata",          # Synonym of Macrophya formosa
        "Tenthredo nimbipennis",       # Valid
        "Lophyrus Abietis",            # Synonym of Neodiprion abietis
        "Lophyrus fulva",              # Synonym of Neodiprion fulviceps
        "Xyela ferruginea",            # Valid
        "Xyela aenea",                 # Invalid (homonym), specimen is Pleroneura aldrichi
        "Tremex columba"               # Valid
    ]

    # A set of the base names that are now considered synonyms for quick lookup.
    synonym_names = {
        "Abia Kennicotti",
        "Ptenos texanus",
        "Ptenos niger",
        "Ptenos nigropectus",
        "Hylotoma abdominalis",
        "Hylotoma miniata",
        "Hylotoma rubiginosa",
        "Nematus chloreus",
        "Emphytus Bollii",
        "Hemichroa fraternalis",
        "Selandria inaequidens",
        "Selandria albicollis",
        "Macrophya excavata",
        "Lophyrus Abietis",
        "Lophyrus fulva",
        "Xyela aenea"
    }

    synonym_indices = []
    for i, full_name in enumerate(species_list_1872):
        # Clean the name to get the base species name for checking
        # e.g., "Cimbex americana, var. Ulmi" -> "Cimbex americana"
        base_name = full_name.split(',')[0].strip()
        
        if base_name in synonym_names:
            # The index is i + 1 because the original list is 1-based.
            synonym_indices.append(i + 1)
            
    # Print the final result as a comma-separated string
    print(",".join(map(str, synonym_indices)))

find_synonym_indices()