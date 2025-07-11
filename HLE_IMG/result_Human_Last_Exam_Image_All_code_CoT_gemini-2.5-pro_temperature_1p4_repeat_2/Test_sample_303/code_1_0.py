def print_poem():
    """
    This function prints the full text of the poem "Ascending Phoenix Terrace in Jinling" by Li Bai.
    """
    poem = {
        "title": "登金陵凤凰台",
        "author": "李白 (唐)",
        "lines": [
            ("凤凰台上凤凰游，凤去台空江自流。", "Fènghuáng tái shàng fènghuáng yóu, fèng qù tái kōng jiāng zì liú.", "On Phoenix Terrace, phoenixes roamed; one left, the terrace empty, the river flows on its own."),
            ("吴宫花草埋幽径，晋代衣冠成古丘。", "Wú gōng huācǎo mái yōujìng, Jìn dài yīguān chéng gǔqiū.", "Wu Palace's flowers and plants lie buried in desolate paths; Jin Dynasty's finery have become ancient mounds."),
            ("三山半落青天外，二水中分白鹭洲。", "Sān shān bàn luò qīngtiān wài, èr shuǐ zhōng fēn Báilù zhōu.", "The Three Mountains fall halfway beyond the blue sky; the two waters are divided by the White Egret Isle."),
            ("总为浮云能蔽日，长安不见使人愁。", "Zǒng wèi fúyún néng bì rì, Cháng'ān bú jiàn shǐ rén chóu.", "It is always the floating clouds that can obscure the sun; not seeing Chang'an makes one feel sorrowful.")
        ]
    }

    print(f"《{poem['title']}》 - {poem['author']}")
    print("-" * 30)
    for line, pinyin, translation in poem['lines']:
        print(line)
        # The following lines can be uncommented to show Pinyin and translation
        # print(pinyin)
        # print(f"({translation})")
        # print()

print_poem()