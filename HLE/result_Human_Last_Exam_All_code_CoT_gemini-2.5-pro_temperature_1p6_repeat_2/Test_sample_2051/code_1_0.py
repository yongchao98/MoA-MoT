#
# This script calculates the degree of the volume polynomial Z_g,n
# and provides the answers to the user's questions.
#

def solve_moduli_questions():
    """
    Answers the two-part question about the volume of the moduli space Z.
    """
    # Part (a): Continuity
    # In the context of Weil-Petersson volumes of moduli spaces, the volume
    # function Z is known to be continuous.
    answer_a = "Yes"

    # Part (b): Degree Calculation
    # We are given the case for g=0, n_+=3, n_-=1.
    g = 0
    n_plus = 3
    n_minus = 1

    # The total number of boundaries 'n' is the sum of n_+ and n_-.
    n = n_plus + n_minus

    # The degree of the polynomial Z_g,n is given by the formula 6g - 6 + 2n.
    degree = 6 * g - 6 + 2 * n

    # Now, we print the answers in a structured format.
    print("(a) Does the property of piecewise polynomiality of Z imply that it is continuous?")
    print(f"Answer: {answer_a}\n")

    print("(b) For g=0, n+=3, and n-=1, what is the degree of the polynomial Z?")
    print(f"The total number of boundaries is n = n_+ + n_- = {n_plus} + {n_minus} = {n}.")
    print(f"The degree of the polynomial is given by the formula: 6g - 6 + 2n.")
    print(f"Substituting g = {g} and n = {n}:")
    print(f"Degree = 6*({g}) - 6 + 2*({n}) = {6*g} - 6 + {2*n} = {degree}")
    print(f"Answer: {degree}")

if __name__ == '__main__':
    solve_moduli_questions()
