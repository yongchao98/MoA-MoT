def solve():
    """
    Identifies and lists the authors of the nine artworks.
    """
    authors = [
        "Han Meilin",
        "Qi Baishi",
        "Luo Zhongli",
        "Mi Fu",
        "Liu Kuo-sung",
        "Zhao Mengfu",
        "Su Shi",
        "Jin Shangyi",
        "Huang Zhou"
    ]

    print("From left to right, top to bottom, the authors are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

    # The final answer in the required format
    final_answer = ", ".join(authors)
    print(f"\n<<<From left to right, top to bottom: {final_answer}>>>")

solve()