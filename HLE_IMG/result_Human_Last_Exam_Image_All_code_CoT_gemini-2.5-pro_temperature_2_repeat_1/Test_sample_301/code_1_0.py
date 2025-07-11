def solve():
    """
    This function identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Liao Bingxiong (廖冰兄)",
        "Qi Baishi (齐白石)",
        "Huang Zhou (黄胄)",
        "Wang Xizhi (王羲之)",
        "Wu Guanzhong (吴冠中)",
        "Li Yong (李邕)",
        "Huaisu (懷素)",
        "He Duoling (何多苓)",
        "Lin Fengmian (林风眠)"
    ]

    print("The authors of the 9 works from left to right, top to bottom are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

solve()
<<<Liao Bingxiong, Qi Baishi, Huang Zhou, Wang Xizhi, Wu Guanzhong, Li Yong, Huaisu, He Duoling, Lin Fengmian>>>