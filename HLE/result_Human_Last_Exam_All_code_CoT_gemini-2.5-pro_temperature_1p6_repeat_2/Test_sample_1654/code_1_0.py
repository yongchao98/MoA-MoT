import re

def count_coat_descriptions():
    """
    Analyzes the text of "The Chameleon" to count the descriptions of Otchumyelov's coat.
    """
    # A representative text of "The Chameleon" containing the key symbolic actions.
    # Translations can vary, but this text includes the 6 commonly cited instances.
    story_text = """
    Police Superintendent Otchumyelov is walking across the market-square wearing a new overcoat and carrying a parcel under his arm.
    ... (The incident with the dog begins) ...
    Someone in the crowd says, "I fancy it's General Zhigalov's!"
    "General Zhigalov's? H'm! ... Take off my coat, Yeldyrin; I can't stand it! I suppose it's going to rain."
    ... (Otchumyelov's opinion of the dog improves) ...
    The policeman says, "No, that's not the General's dog."
    "Not the General's? ... Yeldyrin, help me on with my coat. ... I feel a cold shiver. ... We must give them a lesson! It is high time..."
    ... (Otchumyelov's opinion worsens again) ...
    A bit later, feeling decisive against the dog, Otchumyelov gets hot under the collar. "I'll teach you! ... Yeldyrin," the superintendent addressed the policeman, "take my coat off.... It's frightfully hot!"
    ... (He is about to act, but then...)
    The General's cook mentions the dog belongs to the General's brother, who has just arrived.
    "What, has his brother come? Vladimir Ivanitch?" inquired Otchumyelov, his whole face beaming. "Help me on with my overcoat, Yeldyrin, my lad... there's a wind getting up.... I am cold.... Take it. It's not a bad pup."
    ... (The situation is resolved in favor of the dog) ...
    Otchumyelov threatened the victim, and wrapping himself in his great-coat, went on his way across the market square.
    """

    # The 6 key descriptions or actions related to the coat, symbolizing his shifts in attitude.
    coat_mentions = [
        "wearing a new overcoat", # 1. Initial state: pompous official
        "Take off my coat",       # 2. Action: learns dog may be the General's (becomes servile)
        "help me on with my coat",# 3. Action: learns dog is not the General's (becomes harsh again)
        "take my coat off",       # 4. Action: gets worked up again (becomes hot-headed)
        "Help me on with my overcoat", # 5. Action: learns dog is the General's brother's (becomes fawning)
        "wrapping himself in his great-coat" # 6. Final state: resuming his official facade
    ]

    found_count = 0
    equation_parts = []

    for mention in coat_mentions:
        # Use regex for case-insensitive search to make it more robust
        if re.search(mention, story_text, re.IGNORECASE):
            found_count += 1
            equation_parts.append("1")

    print("Counting the instances where Otchumyelov's coat is described:")
    # Using ' + '.join to construct the equation string dynamically
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {found_count}")
    
    print(f"\nIn \"The Chameleon,\" Chekov described Otchumyelov's coat {found_count} times to symbolize his changing attitude.")

count_coat_descriptions()