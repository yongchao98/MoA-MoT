def find_poem():
    """
    This function prints the complete content of the poem identified from the calligraphy.
    The calligraphy shows the last two lines of Wang Wei's poem "Birdsong Brook".
    """
    
    title = "《鸟鸣涧》"
    author = "【唐】王维"
    line1 = "人闲桂花落，"
    line2 = "夜静春山空。"
    line3 = "月出惊山鸟，"
    line4 = "时鸣春涧中。"
    
    explanation = (
        "这幅书法作品书写的是唐代诗人王维的著名诗作《鸟鸣涧》。\n"
        "图片中的文字是这首诗的后两句：“月出惊山鸟，时鸣春涧中。”\n"
        "这首诗描绘了春山月夜的宁静与生机，意境幽美，是王维“诗中有画”的代表作之一。"
    )

    print("这幅书法作品的完整内容是：")
    print("="*30)
    print(f"{title:^24}")
    print(f"{author:^22}")
    print()
    print(f"{line1:^22}")
    print(f"{line2:^22}")
    print(f"{line3:^22}")
    print(f"{line4:^22}")
    print("="*30)
    print("\n作品介绍：")
    print(explanation)

# Execute the function to print the poem's content
find_poem()