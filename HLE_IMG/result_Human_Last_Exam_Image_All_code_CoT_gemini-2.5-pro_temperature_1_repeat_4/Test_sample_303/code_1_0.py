def display_famous_poem():
    """
    This function prints the full text of the famous Tang dynasty poem "Guan Shan Yue"
    by Li Bai, identified from the calligraphy image.
    """
    title = "《关山月》"
    author = "唐·李白"
    poem_lines = [
        "明月出天山，苍茫云海间。",
        "长风几万里，吹度玉门关。",
        "汉下白登道，胡窥青海湾。",
        "由来征战地，不见有人还。",
        "戍客望边邑，思归多苦颜。",
        "高楼当此夜，叹息未应闲。"
    ]

    print(title)
    print(author)
    print("")  # Add a blank line for spacing
    for line in poem_lines:
        print(line)

# Execute the function to display the poem
display_famous_poem()