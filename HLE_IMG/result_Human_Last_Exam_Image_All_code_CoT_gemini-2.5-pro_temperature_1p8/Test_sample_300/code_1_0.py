def identify_artists():
    """
    This function identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Xu Beihong",       # 1. Top Left: Identified by signature and style.
        "Wu Guanzhong",     # 2. Top Middle: Identified by abstract landscape style.
        "Mao Zedong",       # 3. Top Right: Identified by unique calligraphy style.
        "Qi Baishi",        # 4. Middle Left: Identified by "Baishi" seal and style.
        "Jin Shangyi",      # 5. Middle Middle: A famous oil painting masterpiece.
        "Lin Fengmian",     # 6. Middle Right: Identified by signature artistic style.
        "Wu Guanzhong",     # 7. Bottom Left: Identified by style and signature.
        "Li Keran",         # 8. Bottom Middle: Identified by "Li-style landscape".
        "Qi Baishi"         # 9. Bottom Right: Identified by famous subject matter and style.
    ]

    print("From left to right, top to bottom, the authors of the 9 works are:")
    for i, author in enumerate(authors):
        print(f"Artwork {i+1}: {author}")

identify_artists()