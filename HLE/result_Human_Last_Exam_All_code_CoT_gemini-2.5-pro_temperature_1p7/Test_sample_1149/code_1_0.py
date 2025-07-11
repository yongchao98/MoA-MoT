import sys

def get_lacrimal_recess_depressions():
    """
    Provides the number of depressions in the lacrimal recess of Irritator challengeri.

    This information is based on paleontological findings, specifically the
    description by Sues et al. (2002), which notes the presence of two
    depressions (fossae) within the lacrimal recess. This script presents
    this factual data.
    """
    # The number of depressions as documented in scientific literature.
    num_depressions = 2

    # A factor for the required equation format.
    factor = 1

    # Calculate the result for the equation.
    result = factor * num_depressions

    # Print the answer in a full sentence.
    print(f"The lacrimal recess in Irritator challengeri contains {result} depressions.")

    # Print the equation as requested.
    print(f"Final equation: {factor} * {num_depressions} = {result}")

if __name__ == "__main__":
    get_lacrimal_recess_depressions()
