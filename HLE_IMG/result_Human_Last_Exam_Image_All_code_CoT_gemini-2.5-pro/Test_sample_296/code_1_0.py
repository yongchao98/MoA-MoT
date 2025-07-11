import sys
import io

# Set stdout to utf-8 to handle Cyrillic characters correctly
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def solve_idiom():
    """
    This function explains the reasoning behind identifying the Russian idiom
    and prints the final answer.
    """
    # The first drawing (in red) shows a literal interpretation: "клин вышибают"
    # which means "to knock out a wedge".
    literal_part_1 = "Клин"  # Wedge
    literal_part_2 = "вышибают"  # They knock out

    # The second drawing shows the figurative meaning. A modern taxi (a new "wedge" in the market)
    # is replacing an old horse-drawn cab (an old "wedge"). This illustrates the concept of
    # replacing something with a similar thing.
    figurative_meaning = "A new thing replacing a similar old thing."

    # Combining the literal action with the concept of using a similar object,
    # we get the full idiom.
    idiom_part_1 = "Клин"
    idiom_part_2 = "клином" # "with a wedge" (instrumental case)
    idiom_part_3 = "вышибают"

    print(f"The Russian idiom is formed by combining the parts: '{idiom_part_1}' '{idiom_part_2}' '{idiom_part_3}'.")
    print("Full idiom: Клин клином вышибают")

solve_idiom()