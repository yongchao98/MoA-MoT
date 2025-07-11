def print_poem():
    """
    This function prints the full text of the Tang dynasty poem "The Golden-threaded Robe".
    """
    poem_title = "《金缕衣》"
    author = "（唐）杜秋娘"
    lines = [
        "劝君莫惜金缕衣，",
        "劝君惜取少年时。",
        "花开堪折直须折，",
        "莫待无花空折枝。"
    ]
    pinyin = [
        "quàn jūn mò xī jīn lǚ yī,",
        "quàn jūn xī qǔ shào nián shí.",
        "huā kāi kān zhé zhí xū zhé,",
        "mò dài wú huā kōng zhé zhī."
    ]
    
    print(f"{poem_title} {author}\n")
    for i in range(len(lines)):
        print(f"{lines[i]} ({pinyin[i]})")

print_poem()