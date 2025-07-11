import sys

# This script is designed to run in a simple environment
# where we can print directly to the console.

def solve_word_derivation():
    """
    This function traces a hypothetical word from Proto-Indo-European
    to Middle English, printing each stage of its evolution.
    """
    # Define the forms at each stage of the language evolution.
    pie_form = "*kʷoyséyeti"
    pgmc_form = "*hwaiziþi"
    pwgmc_form = "*hwairiþ"
    oe_form = "hwǣreþ"
    me_form = "whereth"

    print("This script calculates the Middle English reflex of the PIE root *kʷeys-.")
    print("The derivation follows an o-grade causative path.")
    print("-" * 40)
    print("Derivation Steps:")
    print(f"1. Proto-Indo-European: {pie_form}")
    print(f"   (o-grade *kʷoys- + causative suffix -éye- + 3sg ending -ti)")
    print()
    print(f"2. Proto-Germanic:      {pgmc_form}")
    print(f"   (Grimm's & Verner's Laws, vowel shifts)")
    print()
    print(f"3. Proto-West Germanic: {pwgmc_form}")
    print(f"   (Rhotacism of *z > *r)")
    print()
    print(f"4. Old English:         {oe_form}")
    print(f"   (i-umlaut of *ai > *ǣ)")
    print()
    print(f"5. Middle English:      {me_form}")
    print(f"   (Vowel shift *ǣ > ē, spelling change hw- > wh-)")
    print("-" * 40)

    # The prompt asks to "output each number in the final equation!".
    # This is interpreted as showing the full chain of derivation, which looks
    # like a word equation.
    print("The complete 'word equation' is:")
    print(f"{pie_form} > {pgmc_form} > {pwgmc_form} > {oe_form} > {me_form}")
    print()

    print("The final predicted form in Middle English is:")
    # The final answer is printed clearly.
    final_answer = me_form
    print(final_answer)

solve_word_derivation()
# The final answer is also enclosed in the special markers as requested.
sys.stdout.write("<<<whereth>>>")