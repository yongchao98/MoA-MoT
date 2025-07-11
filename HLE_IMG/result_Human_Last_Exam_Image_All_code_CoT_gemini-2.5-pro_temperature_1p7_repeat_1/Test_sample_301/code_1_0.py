def solve():
    """
    This function identifies and lists the authors of the nine artworks.
    """
    authors = [
        "Han Meilin",
        "Qi Baishi",
        "Ai Xuan",
        "Huaisu",
        "Wu Guanzhong",
        "Chu Suiliang",
        "Zhang Xu",
        "Jin Shangyi",
        "Lin Fengmian"
    ]

    print("The authors of the nine works from left to right, top to bottom are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

    final_answer = ", ".join(authors)
    print(f"\n<<<1. Han Meilin, 2. Qi Baishi, 3. Ai Xuan, 4. Huaisu, 5. Wu Guanzhong, 6. Chu Suiliang, 7. Zhang Xu, 8. Jin Shangyi, 9. Lin Fengmian>>>")

solve()