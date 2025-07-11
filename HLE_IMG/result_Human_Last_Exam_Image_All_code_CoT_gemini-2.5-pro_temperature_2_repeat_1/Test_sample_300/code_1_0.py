def solve():
    """
    This function lists the authors of the nine artworks.
    """
    authors = [
        # Top Row
        "Huang Zhou (黄胄)", 
        "Wu Guanzhong (吴冠中)", 
        "Huaisu (怀素)",
        
        # Middle Row
        "Mao Zedong (毛泽东)", 
        "Jin Shangyi (靳尚谊)", 
        "Lin Fengmian (林风眠)",
        
        # Bottom Row
        "Wu Guanzhong (吴冠中)", 
        "Li Keran (李可染)", 
        "Qi Baishi (齐白石)"
    ]
    
    print("From left to right, top to bottom, the authors are:")
    for i, author in enumerate(authors):
        print(f"{i+1}. {author}")

solve()