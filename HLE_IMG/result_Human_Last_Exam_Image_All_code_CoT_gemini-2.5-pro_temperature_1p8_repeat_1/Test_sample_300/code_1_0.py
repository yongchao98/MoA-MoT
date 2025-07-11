def solve():
    """
    Identifies and lists the authors of the 9 artworks.
    """
    authors = [
        "Huang Zhou (黄胄)",
        "Pan Tianshou (潘天寿)",
        "Xu Wei (徐渭)",
        "Fu Shan (傅山)",
        "Chen Danqing (陈丹青)",
        "Lin Fengmian (林風眠)",
        "Wu Guanzhong (吴冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齐白石)"
    ]

    print("From left to right, from top to bottom, the authors are:")
    for i, author in enumerate(authors, 1):
        print(f"Work {i}: {author}")

solve()