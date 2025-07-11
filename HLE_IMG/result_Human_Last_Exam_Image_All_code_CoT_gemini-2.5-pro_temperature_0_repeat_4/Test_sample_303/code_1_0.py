def print_poem():
    """
    This function prints the full text of the poem "Spring View" (春望) by Du Fu.
    """
    title = "春望"
    author = "杜甫 (唐)"
    
    poem = [
        "国破山河在，",
        "城春草木深。",
        "感时花溅泪，",
        "恨别鸟惊心。",
        "烽火连三月，",
        "家书抵万金。",
        "白头搔更短，",
        "浑欲不胜簪。"
    ]
    
    print(f"作品名称：《{title}》")
    print(f"作者：{author}")
    print("--------------------")
    for line in poem:
        print(line)
    print("--------------------")

print("这幅书法作品的完整内容是唐代诗人杜甫的《春望》：\n")
print_poem()