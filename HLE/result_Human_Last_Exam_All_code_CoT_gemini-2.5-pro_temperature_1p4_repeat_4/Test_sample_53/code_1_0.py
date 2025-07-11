# This script derives the hypothetical Middle English word step-by-step.

# Step 1: Proto-Indo-European (PIE)
pie_root = "*kʷeys-"
pie_causative = "*kʷoyséye-"
print(f"1. PIE Root & Causative")
print(f"   - The root is {pie_root} ('to see, to heed').")
print(f"   - The o-grade causative form is {pie_causative} ('to show').")
print("-" * 40)

# Step 2: Proto-Germanic (PGmc)
pgmc_form = "*hʷaizjaną"
print(f"2. Development into Proto-Germanic")
print(f"   - Grimm's Law (*kʷ > *hʷ) + Verner's Law (*s > *z) + Vowel/Suffix changes.")
print(f"   - The resulting infinitive is {pgmc_form}.")
print("-" * 40)

# Step 3: Old English (OE)
oe_3sg_present = "hǣseþ"
print(f"3. Development into Old English")
print(f"   - i-Umlaut (*ai > *ǣ) + Consonant shifts (*hʷ > *h, *zj > s).")
print(f"   - Adding the 3rd person singular present ending '-eþ'.")
print(f"   - The resulting OE form is: {oe_3sg_present}.")
print("-" * 40)

# Step 4: Middle English (ME)
me_3sg_present = "heseth"
print(f"4. Development into Middle English")
print(f"   - Vowel shift (OE *ǣ > ME long open ē).")
print(f"   - Ending preservation (OE -eþ > ME -eth).")
print(f"   - The final ME form meaning 'he shows' is:")
print(f"   -> {me_3sg_present}")
print("-" * 40)

# Final Answer
print(f"<<<{me_3sg_present}>>>")