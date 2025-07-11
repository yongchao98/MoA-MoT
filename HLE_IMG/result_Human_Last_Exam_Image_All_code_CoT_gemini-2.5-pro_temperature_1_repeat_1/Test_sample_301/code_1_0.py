def solve():
    """
    This function identifies and lists the authors of the 9 artworks.
    """
    authors = [
        "Han Meilin (韩美林)",
        "Qi Baishi (齐白石)",
        "Ai Xuan (艾轩)",
        "Huang Tingjian (黄庭坚)",
        "Wu Guanzhong (吴冠中)",
        "Chu Suiliang (褚遂良)",
        "Huaisu (怀素)",
        "Chen Yifei (陈逸飞)",
        "He Jiaying (何家英)"
    ]

    print("The authors of the 9 artworks, from left to right and top to bottom, are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

solve()