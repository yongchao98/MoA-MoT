import math

def solve_overhang_puzzle():
    """
    Solves the 3-cube maximal overhang problem.

    The problem asks for integers a, b, c for the maximal overhang expression
    (a + sqrt(b)) / (1 + c).

    The maximal overhang for three identical cubes is a known result from computational
    geometry and physics, achieved by a configuration involving rotation.
    One cube is placed on the table, rotated by an optimal angle (30 degrees).
    The other two cubes are placed on top of this rotated base.

    The value of this maximal overhang is (1 + sqrt(3)) / 2 side-lengths.

    We need to match this value to the format (a + sqrt(b)) / (1 + c),
    with the given constraints for integers a, b, and c.

    - a, b, c must be non-negative integers.
    - sqrt(b) must be either zero or non-integer.
    - c must be minimal.

    Let's compare the maximal overhang with the given format:
    Maximal Overhang = (1 + sqrt(3)) / 2

    Given Format = (a + sqrt(b)) / (1 + c)

    By comparing the numerators:
    a + sqrt(b) = 1 + sqrt(3)
    This gives a = 1 and b = 3.

    By comparing the denominators:
    1 + c = 2
    This gives c = 1.

    Let's check the constraints for a=1, b=3, c=1:
    - a, b, c are non-negative integers: 1, 3, 1 are all valid.
    - sqrt(b) = sqrt(3) is non-integer, which is valid.
    - c = 1 is the minimal non-negative integer for the denominator to be 2.

    So, the solution is a=1, b=3, c=1.

    The code will now print the final equation and the answer.
    """
    a = 1
    b = 3
    c = 1

    print("The maximal overhang for three cubes can be expressed as (a + sqrt(b)) / (1 + c).")
    print(f"Based on the known optimal configuration, the value is ({a} + sqrt({b})) / ({1} + {c}).")
    
    overhang_val = (a + math.sqrt(b)) / (1 + c)
    print(f"This evaluates to approximately {overhang_val:.4f} side-lengths.")
    
    # Final answer in the specified format
    final_answer = f"{a} {b} {c}"
    print("\nThe values for a, b, and c are:")
    print(final_answer)

solve_overhang_puzzle()
# <<<1 3 1>>>