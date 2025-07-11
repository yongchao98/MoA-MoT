def reveal_poem():
    """
    This function prints the identified poem in Simplified Chinese.
    The poem is read from the calligraphy image, which is in bird-and-worm seal script.
    """
    poem_line_1 = "闺中思"  # English: Thoughts from the boudoir
    poem_line_2 = "君出征"  # English: My lord goes on a campaign
    poem_line_3 = "飞霜龙"  # English: (with his) Flying Frost Dragon (sword)
    
    full_poem = f"{poem_line_1}，{poem_line_2}，{poem_line_3}。"
    
    print("这幅作品的完整内容是：")
    print(full_poem)
    print("\n大致释义为：")
    print("在闺房中深深地思念，我的夫君已出征远行，（伴随他的是名为）飞霜龙的宝剑。")

reveal_poem()