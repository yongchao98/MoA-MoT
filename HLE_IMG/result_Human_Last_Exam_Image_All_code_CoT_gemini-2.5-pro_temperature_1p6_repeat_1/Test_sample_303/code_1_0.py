# -*- coding: utf-8 -*-

def print_poem():
    """
    This function prints the full text of the poem "Entering Ruoye Stream" by Wang Ji.
    """
    poem_title = "《入若耶溪》"
    author_info = "  （南朝·梁）王籍"
    line1 = "艅艎何泛泛，空水共悠悠。"
    line2 = "阴霞生远岫，阳景逐回流。"
    line3 = "蝉噪林愈静，鸟鸣山更幽。"
    line4 = "此地动归念，长年悲倦游。"
    
    full_poem = f"""
这幅书法作品节选了南朝诗人王籍的名篇《入若耶溪》。
作品中的两句是千古名句：“蝉噪林愈静，鸟鸣山更幽”。

以下是这首诗的全文（简体中文）：

{poem_title}
{author_info}

{line1}
{line2}
{line3}
{line4}
"""
    print(full_poem)

print_poem()