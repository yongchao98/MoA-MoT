def find_evaluation():
    """
    This function retrieves and displays Liu Ban's evaluation of Zhang Ji's seven-character poems.
    """
    # Key figures and topic
    evaluator = "刘攽"
    poet = "张籍"
    poet_courtesy_name = "张文昌"  # 张籍的字是文昌，古人常以字相称
    poem_style = "七言诗"

    # The exact evaluation found in classical literary criticism texts.
    # This famous quote is widely attributed to Liu Ban.
    # The evaluation is broken down into its key components for clarity.
    part_1_praise = "多佳句"
    part_2_criticism = "殊未及格"
    part_3_reason = "格律不精耳"

    # Construct the full original quote
    full_quote = f"“{poet_courtesy_name}{poem_style}，{part_1_praise}，{part_2_criticism}，所以然者，{part_3_reason}。”"

    # Print the final answer clearly
    print(f"{evaluator}对{poet}的{poem_style}的精准评价原文是：")
    print(full_quote)

    # Provide a simple explanation in modern Chinese
    print("\n---")
    print("释义：")
    print(f"这段评价的意思是：张籍（字文昌）的七言诗，虽然有很多（'{part_1_praise}'）的好句子，")
    print(f"但是整体上还远远不够标准（'{part_2_criticism}'）。")
    print(f"造成这种情况的原因，是其诗歌的声韵格律不够精通和严谨（'{part_3_reason}'）。")

if __name__ == "__main__":
    find_evaluation()