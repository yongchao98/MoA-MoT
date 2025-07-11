def find_authors():
    """
    Identifies and prints the authors of the 9 artworks in the image grid.
    The order is from left to right, top to bottom.
    """
    authors = [
        "1. (Top-Left): Luo Ping (罗平)",
        "2. (Top-Middle): Qi Baishi (齐白石)",
        "3. (Top-Right): Jin Shangyi (靳尚谊)",
        "4. (Middle-Left): Mi Fu (米芾)",
        "5. (Middle-Middle): Fong Chung-Ray (馮鍾睿)",
        "6. (Middle-Right): Chu Suiliang (褚遂良)",
        "7. (Bottom-Left): Mi Fu (米芾)",
        "8. (Bottom-Middle): Jin Shangyi (靳尚谊)",
        "9. (Bottom-Right): He Jiaying (何家英)"
    ]

    print("The authors of the nine works are:")
    for author in authors:
        print(author)

find_authors()