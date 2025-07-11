def find_authors():
    """
    This function identifies and prints the authors of the nine artworks
    in the image grid, from left to right, top to bottom.
    """
    authors = {
        1: "Lao Shu (老树)",
        2: "Qi Baishi (齐白石)",
        3: "Wang Zhaomin (王兆民)",
        4: "Mao Zedong (毛泽东)",
        5: "Liu Guosong (刘国松)",
        6: "Wang Xizhi (王羲之)",
        7: "Huaisu (怀素)",
        8: "Jin Shangyi (靳尚谊)",
        9: "Wang Meifang (王美芳)"
    }

    print("The authors of the 9 works are:")
    for i in range(1, 10):
        print(f"{i}. {authors[i]}")

find_authors()