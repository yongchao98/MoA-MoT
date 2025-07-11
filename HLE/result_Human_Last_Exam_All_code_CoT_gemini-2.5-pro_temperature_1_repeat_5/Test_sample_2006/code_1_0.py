import math
from fractions import Fraction

def solve_haiku_riddle():
    """
    Solves the riddle presented in the haiku by treating it as a
    quadratic equation, where coefficients are derived from letter counts.
    """
    # The lines of the haiku
    line1 = "An August tempest"
    line2 = "Twice fifteen brings winds of change"
    line3 = "A divine one yields"

    # Calculate coefficients by summing the length of the words in each line
    a = sum(len(word) for word in line1.split())
    b = sum(len(word) for word in line2.split())
    c = sum(len(word) for word in line3.split())

    print("Based on the riddle, we form the quadratic equation ax^2 + bx + c = 0.")
    print(f"The final equation is: {a}x^2 + {b}x + {c} = 0")
    print("-" * 30)

    # Calculate the discriminant (b^2 - 4ac)
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("The equation has no real roots.")
    else:
        # Use integer square root for precision since the discriminant is a perfect square (1)
        sqrt_discriminant = math.isqrt(discriminant)

        # Calculate the numerators and the common denominator
        numerator1 = -b + sqrt_discriminant
        numerator2 = -b - sqrt_discriminant
        denominator = 2 * a

        # Represent roots as fractions for exactness
        root1 = Fraction(numerator1, denominator)
        root2 = Fraction(numerator2, denominator)

        # Sort the roots in numerical ("alphabetical") order
        answers = sorted([root1, root2])

        print(f"The answers (the 'Bays'), in alphabetical order, are:")
        print(f"Root 1: {answers[0]}")
        print(f"Root 2: {answers[1]}")

solve_haiku_riddle()