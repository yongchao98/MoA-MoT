def solve():
    """
    Identifies the author of each of the 9 artworks from left to right, top to bottom.
    """
    authors = [
        "1. Han Meilin",
        "2. Qi Baishi",
        "3. Jin Shangyi",
        "4. Wang Xizhi",
        "5. Liu Kuo-sung",
        "6. Yan Zhenqing",
        "7. Su Shi",
        "8. Yang Feiyun",
        "9. Lin Yong"
    ]
    
    for author in authors:
        print(author)

    final_answer = ", ".join([name.split('. ')[1] for name in authors])
    # The final answer format is specified as <<<answer content>>>
    # The content should be the list of authors.
    print(f"\n<<<{final_answer}>>>")

solve()