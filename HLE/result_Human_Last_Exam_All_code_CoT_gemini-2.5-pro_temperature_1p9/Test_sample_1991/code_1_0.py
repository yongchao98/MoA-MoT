import sys
# This script will generate and print a phrase that meets the user's specific constraints.

def solve_puzzle():
    """
    Solves the creative puzzle by crafting a phrase according to the rules.

    Rules:
    1. Fourteen syllables.
    2. Describes abnormal cell death in Luys' body causing cancer.
    3. Rhyme and scansion match a Disney song title ("Supercalifragilisticexpialidocious").
    4. Each word capitalized.
    5. No hyphens.
    """

    # The crafted phrase:
    # Syllable count: Luys'(1) + Bo-dy(2) + Le-sion(2) + Pro-lif-er-a-tion(5) + Is(1) + Per-ni-cious(3) = 14 syllables.
    # The phrase describes a cancerous lesion's growth (proliferation) in the Luys' Body.
    # "Pernicious" rhymes with the "-docious" part of the reference song title and matches its rhythm.
    phrase = "Luys' Body Lesion Proliferation Is Pernicious"

    # Capitalize each word in the phrase.
    capitalized_phrase = phrase.upper()

    # The prompt also states: "Remember in the final code you still need to output each number in the final equation!"
    # As there is no equation, this instruction seems non-applicable.
    # The only number in the prompt is "fourteen". I will print the phrase as the final result.

    print(capitalized_phrase)

# Execute the function to print the final answer.
solve_puzzle()