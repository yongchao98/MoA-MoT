import math
import fractions

def solve_haiku_riddle():
    """
    This function solves the riddle presented in the haiku by treating it
    as a quadratic equation.
    """
    # Step 1: Decode the coefficients a, b, and c from the haiku.
    # a: "An August tempest" -> "August" is the 8th month.
    a = 8
    # b: "Twice fifteen" -> 2 * 15.
    b = 30
    # c: "A divine one yields" -> 7 is often considered a "divine" number.
    # This choice makes the discriminant a perfect square, a common feature in such riddles.
    c = 7

    # The prompt requires outputting each number in the final equation.
    print(f"The riddle describes the quadratic equation: {a}x^2 + {b}x + {c} = 0")
    print("-" * 40)

    # Step 2: Solve the quadratic equation using the quadratic formula.
    # x = [-b ± sqrt(b^2 - 4ac)] / 2a
    discriminant = (b**2) - (4 * a * c)

    if discriminant < 0:
        print("The equation has no real roots.")
        return

    # Calculate the two roots
    sqrt_discriminant = math.sqrt(discriminant)
    
    # Root 1
    root1_numerator = -b + sqrt_discriminant
    root1_denominator = 2 * a
    root1 = root1_numerator / root1_denominator
    
    # Root 2
    root2_numerator = -b - sqrt_discriminant
    root2_denominator = 2 * a
    root2 = root2_numerator / root2_denominator

    # Represent roots as precise fractions
    frac1 = fractions.Fraction(root1).limit_denominator()
    frac2 = fractions.Fraction(root2).limit_denominator()
    
    print(f"The two solutions (roots) are {root1} and {root2}.")
    print(f"Expressed as fractions, the roots are {frac1} and {frac2}.")
    print("-" * 40)
    
    # Step 3: Determine the final answer based on the alphabetical order constraint.
    # The word representations of the fractions need to be compared.
    # We assign the word forms based on which fraction corresponds to which root.
    root_words = {
        str(frac1): "negative one quarter",
        str(frac2): "negative seven halves"
    }

    print("The question asks for the answer in alphabetical order.")
    print("Comparing the English names of the roots:")
    print(f"'{frac1}' is read as '{root_words[str(frac1)]}'")
    print(f"'{frac2}' is read as '{root_words[str(frac2)]}'")

    # Determine which comes first alphabetically
    final_answer_word = min(root_words.values())
    
    final_answer_fraction = ""
    for k, v in root_words.items():
        if v == final_answer_word:
            final_answer_fraction = k
            break
            
    print(f"\nAlphabetically, '{final_answer_word}' comes first.")
    print(f"Therefore, the final answer is {final_answer_fraction}.")

solve_haiku_riddle()