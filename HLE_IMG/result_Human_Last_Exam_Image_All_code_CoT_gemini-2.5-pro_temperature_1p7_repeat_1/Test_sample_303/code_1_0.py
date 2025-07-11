def explain_calligraphy():
    """
    Identifies and explains the content of the provided calligraphy image.
    """
    # The source of the text
    source = "《周易·乾·文言》"
    
    # The 9 characters written in the calligraphy
    text_in_work = "云从龙，风从虎，圣人作"
    
    # The full sentence from which the text is excerpted
    full_quote = "云从龙，风从虎，圣人作而万物睹"
    
    # Explanation of the meaning
    meaning = ("这是一个比喻，意思是云跟随着龙而出现，风跟随着虎而刮起，"
               "圣人兴起则天下万物都会彰显其风采。它揭示了同类事物之间会互相感应的道理。")

    print(f"这幅书法作品是用一种叫做“鸟虫篆”的高度艺术化字体书写的。")
    print(f"其内容出自中国古典文献{source}。")
    print(f"作品中的九个字是：“{text_in_work}”。")
    print(f"它所在的完整句子是：“{full_quote}”。")
    print(f"这句话的含义是：{meaning}")

explain_calligraphy()