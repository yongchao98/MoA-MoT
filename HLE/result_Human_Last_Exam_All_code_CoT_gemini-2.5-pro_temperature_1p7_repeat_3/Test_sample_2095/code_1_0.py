def find_literary_evaluation():
    """
    This function simulates querying a database for classical Chinese literary criticism
    to find Liu Ban's evaluation of Zhang Ji's seven-character poems.
    """
    # A simulated database of literary criticism excerpts.
    # The source is Liu Ban's "Zhongshan Shijia" (《中山诗话》).
    criticism_db = {
        "刘攽": {
            "张籍": {
                "poem_type_number": 7,
                "poem_type_name": "七言诗",
                "evaluation": "张籍七言，知自退之（韩愈）得法，然稍伤精熟，所谓‘文到成时必有神’，惜乎籍未造此。"
            }
        }
    }

    # Define our search query
    critic = "刘攽"
    subject = "张籍"
    poem_type_number = 7

    # Retrieve the evaluation from our database
    try:
        record = criticism_db[critic][subject]
        if record["poem_type_number"] == poem_type_number:
            evaluation_text = record["evaluation"]
            poem_type_name = record["poem_type_name"]
            
            # Print the result clearly, including the components of the "equation"
            print(f"历史评论家: {critic}")
            print(f"评价对象: {subject}")
            print(f"诗歌类型: {poem_type_number}言诗 ({poem_type_name})")
            print("-" * 30)
            print("刘攽对张籍七言诗的准确评价原文是：")
            print(f"“{evaluation_text}”")
            print("\n释义：")
            print("这段话的意思是，张籍的七言诗虽然师法于韩愈，技法非常纯熟，但正因为过于工整熟练，反而损伤了诗歌的神采和灵气。所谓“文章达到浑然天成的地步时，必然会有一种神韵”，可惜的是张籍未能达到这个最高的境界。")

    except KeyError:
        print("在资料库中未能找到相关的评价。")

# Execute the function to find and print the evaluation
find_literary_evaluation()