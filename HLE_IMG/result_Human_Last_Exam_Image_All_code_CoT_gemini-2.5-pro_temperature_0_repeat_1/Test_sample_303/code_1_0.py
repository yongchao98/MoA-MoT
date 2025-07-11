def print_poem():
    """
    Prints the full text of the poem "Spring View" (春望) by Du Fu.
    """
    title = "春望"
    author = "杜甫 (唐)"
    poem_lines = [
        "国破山河在，城春草木深。",
        "感时花溅泪，恨别鸟惊心。",
        "烽火连三月，家书抵万金。",
        "白头搔更短，浑欲不胜簪。"
    ]
    
    print(f"《{title}》")
    print(f"{author}\n")
    for line in poem_lines:
        print(line)

print_poem()