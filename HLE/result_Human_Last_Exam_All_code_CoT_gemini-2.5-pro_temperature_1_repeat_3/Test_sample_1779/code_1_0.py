import math

def main():
    """
    This script analyzes the decidability of a language L where a string w is in L
    if and only if its length is a perfect square. It implements a decider algorithm
    to prove that the language is, in fact, decidable.
    """

    print("The problem asks if there exists a decidable language L for a Turing machine T,")
    print("where T halts on an input string w if and only if len(w) is a perfect square.")
    print("\nThis means L is the set of all strings whose length is a perfect square.")
    print("A language is decidable if an algorithm exists that can always determine membership and always halts.")
    print("We can prove L is decidable by constructing such an algorithm.\n")

    def decide_language_L(w: str):
        """
        This function is an algorithm (a decider) for the language L.
        It takes a string 'w', determines if its length is a perfect square,
        and prints the step-by-step logic. It always halts.
        """
        n = len(w)
        print(f"--- Analyzing string: '{w}' ---")
        print(f"Step 1: Find the length of the string.")
        print(f"Length n = {n}")

        print("\nStep 2: Check if n is a perfect square.")
        # An efficient way to check for a perfect square.
        # Note: We handle n=0 as a special case for clarity, though math.isqrt(0) works.
        if n < 0:
            # This case is not possible for string length but is included for completeness.
            is_perfect_square = False
            equation_str = "N/A as length is negative"
        else:
            s = int(math.sqrt(n))
            is_perfect_square = (s * s == n)
            # This part fulfills the requirement to output each number in the final equation.
            equation_str = f"{s} * {s} == {n}"
        
        print(f"To do this, we calculate the integer square root of {n}, which is {s}.")
        print(f"We then test the equation: {equation_str}")
        print(f"The result of this check is: {is_perfect_square}")

        print("\nStep 3: Conclude based on the result.")
        if is_perfect_square:
            print(f"Conclusion: The length {n} is a perfect square. The string '{w}' is in language L.")
        else:
            print(f"Conclusion: The length {n} is NOT a perfect square. The string '{w}' is NOT in language L.")
        print("-" * (len(w) + 24))


    # Demonstrate with a string whose length is a perfect square (9).
    w1 = "111000111"
    decide_language_L(w1)

    print("\n")

    # Demonstrate with a string whose length is not a perfect square (10).
    w2 = "helloworld"
    decide_language_L(w2)

    print("\nBecause we have an algorithm that correctly decides membership in L for any input string")
    print("and is guaranteed to halt, the language L is decidable.")
    print("\nTherefore, the answer to the question is Yes.")

if __name__ == "__main__":
    main()

<<<Yes>>>