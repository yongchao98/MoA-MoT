def derive_middle_english_form():
    """
    This function traces a hypothetical word from its Proto-Indo-European root
    to its Middle English form and prints the derivation.
    """
    # Define the forms at each stage of development
    pie_form = "*kʷoys-eye-"
    pgmc_form = "*hʷais-ija-"
    oe_form = "hwǣseþ"
    me_form = "wheseth"

    # Explain the result
    print("The step-by-step derivation of the hypothetical word is as follows:")

    # Print each form in the derivation path
    print(f"1. Proto-Indo-European (o-grade causative stem): {pie_form}")
    print(f"2. Proto-Germanic (Class I weak verb stem): {pgmc_form}")
    print(f"3. Old English (3rd person sg. present): {oe_form}")
    print(f"4. Middle English (3rd person sg. present): {me_form}")

    print("\nThe final derivation can be summarized in the following equation, showing each form:")
    # The final "equation" with each form as requested.
    print(f"{pie_form} > {pgmc_form} > {oe_form} > {me_form}")

    print("\nThe predicted Middle English form for 'he shows' is:")
    print(me_form)

if __name__ == "__main__":
    derive_middle_english_form()
