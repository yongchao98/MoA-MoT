import textwrap

def solve_mirror_puzzle():
    """
    Explains the logic for solving the Mirror and Oni puzzle and prints the final answer.
    """

    explanation = """
    The best way to expose the demon is to test the fundamental geometry of a reflection, which is something a simple mimic would likely fail to replicate. Option G provides the most elegant and foolproof method.

    Here is the step-by-step process:

    1.  Place a large, opaque object (like a folding screen) on the floor between you and the mirror. This screen should completely block your view of the mirror and the mirror's view of you.

    2.  Take a single step to your right, so that you are peeking out from behind the right edge of the screen.

    3.  Observe the 'reflection'.

    -   Case 1: True Reflection
        A mirror creates a left-right inverted image. When you step to your right, your reflection effectively steps to its left. From your vantage point, you will see your reflection also peeking from the right side of the reflected screen.

    -   Case 2: Demon Mimic
        A demon in a 'mirror world' simply copying your actions would see you step to your right and would, in turn, step to its own right. From your vantage point, the demon would appear on the opposite (left) side of the room, peeking from the left side of the screen.

    This simple positional test creates a clear paradox that the demon cannot resolve without knowing it must simulate reflection physics rather than just copy your actions. This reveals the illusion.
    """

    print(textwrap.dedent(explanation).strip())

solve_mirror_puzzle()
print("<<<G>>>")