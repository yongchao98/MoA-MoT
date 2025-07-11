def solve():
    """
    Identifies the author of each of the 9 works and prints the names.
    """
    authors = [
        "Huang Zhou",
        "Wu Guanzhong",
        "Mao Zedong",
        "Shitao",
        "Jin Shangyi",
        "Lin Fengmian",
        "Wu Guanzhong",
        "Li Keran",
        "Qi Baishi"
    ]
    
    print("From left to right, top to bottom, the authors are:")
    for i, author in enumerate(authors, 1):
        print(f"Work {i}: {author}")
        
    final_answer = ", ".join(authors)
    # The final answer format is specified as <<<answer content>>>
    # The user request asks for a list of authors.
    # So I will present the final answer as a comma-separated string.
    # print(f"\n<<< {final_answer} >>>")

solve()