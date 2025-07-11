def solve():
    """
    This function identifies and prints the authors of the 9 artworks.
    The order is from left to right, top to bottom.
    """
    authors = [
        "Han Meilin (韩美林)",
        "Qi Baishi (齐白石)",
        "Jin Shangyi (靳尚谊)",
        "Su Shi (苏轼)",
        "Liu Guosong (刘国松)",
        "Mi Fu (米芾)",
        "Huaisu (怀素)",
        "Yang Feiyun (杨飞云)",
        "Lin Fengmian (林风眠)"
    ]

    print("The authors of the artworks from left to right, top to bottom are:")
    for i, author in enumerate(authors, 1):
        print(f"Work {i}: {author}")

solve()
<<<
Han Meilin (韩美林)
Qi Baishi (齐白石)
Jin Shangyi (靳尚谊)
Su Shi (苏轼)
Liu Guosong (刘国松)
Mi Fu (米芾)
Huaisu (怀素)
Yang Feiyun (杨飞云)
Lin Fengmian (林风眠)
>>>