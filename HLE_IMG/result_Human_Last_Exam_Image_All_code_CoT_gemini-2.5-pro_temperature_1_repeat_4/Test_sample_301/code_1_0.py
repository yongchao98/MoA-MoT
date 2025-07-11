def solve_art_authors():
    """
    Identifies and prints the author of each of the 9 artworks in the grid,
    from left to right, top to bottom.
    """
    authors = {
        1: "Lao Shu (老树)",
        2: "Qi Baishi (齐白石)",
        3: "Huang Zhongyang (黄中羊)",
        4: "Huang Tingjian (黄庭坚)",
        5: "Wu Guanzhong (吴冠中)",
        6: "Chu Suiliang (褚遂良)",
        7: "Huai Su (怀素)",
        8: "Jin Shangyi (靳尚谊)",
        9: "He Jiaying (何家英)"
    }

    print("From left to right, from top to bottom, the authors of the 9 works are:")
    for i in range(1, 10):
        print(f"{i}. {authors[i]}")

solve_art_authors()