def solve():
    """
    This function identifies and lists the authors of the nine artworks.
    """
    authors = [
        "Han Meilin (韩美林)",
        "Qi Baishi (齐白石)",
        "Huang Zhou (黄胄)",
        "Huaisu (怀素)",
        "Wu Guanzhong (吴冠中)",
        "Chu Suiliang (褚遂良)",
        "Wang Xizhi (王羲之)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)"
    ]

    print("From left to right, top to bottom, the authors are:")
    for i, author in enumerate(authors):
        print(f"{i+1}. {author}")

    # Final answer format
    final_answer = ", ".join(authors)
    print(f"\n<<<From left to right, top to bottom: {final_answer}>>>")

solve()