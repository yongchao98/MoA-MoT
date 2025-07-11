def find_authors():
    """
    Identifies and prints the authors of the 9 artworks in the image.
    The order is from left to right, top to bottom.
    """
    authors = [
        "1. Lao Shu (老树)",
        "2. Qi Baishi (齐白石)",
        "3. Ai Xuan (艾轩)",
        "4. Huai Su (懷素)",
        "5. Walasse Ting (丁雄泉)",
        "6. Wang Xizhi (王羲之)",
        "7. Huai Su (懷素)",
        "8. Jin Shangyi (靳尚谊)",
        "9. Lin Fengmian (林风眠)"
    ]

    print("The authors of the 9 works, from left to right, top to bottom are:")
    for author in authors:
        print(author)

find_authors()