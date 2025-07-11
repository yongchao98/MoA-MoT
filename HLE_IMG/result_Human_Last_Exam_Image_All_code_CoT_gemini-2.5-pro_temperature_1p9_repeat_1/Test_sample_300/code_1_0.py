def solve():
    """
    This function identifies and lists the artists of the nine artworks.
    """
    artists = [
        "Huang Zhou (黄胄)", "Li Keran (李可染)", "Mao Zedong (毛泽东)",
        "Qi Gong (启功)", "Jin Shangyi (靳尚谊)", "Lin Fengmian (林风眠)",
        "Wu Guanzhong (吴冠中)", "Fu Baoshi (傅抱石)", "Qi Baishi (齐白石)"
    ]

    print("From left to right, from top to bottom, the authors are:")
    for i, artist in enumerate(artists, 1):
        print(f"Work {i}: {artist}")

solve()