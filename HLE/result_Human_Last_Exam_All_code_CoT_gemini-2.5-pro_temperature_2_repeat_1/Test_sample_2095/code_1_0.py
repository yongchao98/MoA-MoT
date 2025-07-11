def liu_ban_evaluation():
    """
    This function prints Liu Ban's evaluation of Zhang Ji's seven-character poems.
    """
    # The evaluation comes from Liu Ban's "Zhongshan Shihua" (《中山诗话》).
    evaluation_source = "北宋刘攽在其《中山诗话》中对张籍的七言诗有如下评价："

    # The original quote from the text.
    original_quote = "“张籍七言诗，拙涩少味，独《征西将》一篇，稍可与《燕歌行》相上下耳。”"

    # Explanation of the evaluation.
    explanation = (
        "这段评价的意思是：\n"
        "刘攽认为张籍的七言诗总体上写得比较笨拙、晦涩，缺乏韵味。\n"
        "但是，他特别指出了一个例外，就是《征西将》这首诗。他认为这唯一的一篇佳作，"
        "其水平几乎可以和高适的名篇《燕歌行》不相上下了。"
    )

    # Combine all parts for the final answer.
    full_answer = f"{evaluation_source}\n{original_quote}\n\n{explanation}"

    print(full_answer)

# Execute the function to display the answer.
liu_ban_evaluation()

final_answer = "刘攽在《中山诗话》中评价张籍的七言诗为：“张籍七言诗，拙涩少味，独《征西将》一篇，稍可与《燕歌行》相上下耳。” 意思是，他认为张籍的七言诗总体上写得比较笨拙晦涩、缺乏味道，只有一个例外，就是《征西将》这首诗写得非常好，水平几乎可以和高适的名篇《燕歌行》相媲美。"
print(f"\n<<<{final_answer}>>>")