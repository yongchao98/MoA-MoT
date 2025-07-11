import re

def find_coat_descriptions():
    """
    This function finds and counts the descriptions of Otchumyelov's
    coat actions in Anton Chekov's "The Chameleon".
    """
    # Full text of "The Chameleon" by Anton Chekov (trans. Constance Garnett)
    story_text = """
    Police Superintendent Otchumyelov is walking across the market square wearing a new overcoat and carrying a parcel under his arm. A red-haired policeman strides after him with a sieve full of confiscated gooseberries in his hands. There is silence all around. Not a soul in the square.... The open doors of the shops and taverns look out upon God's world disconsolately, like hungry mouths; there is not even a beggar near them.

    "So you bite, you damned brute?" Otchumyelov hears suddenly. "Lads, don't let him go! Biting is prohibited nowadays! Hold him! Ah!... ah!"

    There is the sound of a dog yelping. Otchumyelov looks in the direction of the sound and sees a dog, hopping on three legs and looking about her, run out of Pitchugin's timber-yard. A man in a starched cotton shirt, with his waistcoat unbuttoned, is chasing her. He runs after her, and throwing his body forward falls down and seizes the dog by her hind legs. The sound of yelping and the cry of "Don't let go!" are heard once more. Sleepy countenances are protruded from the shops, and soon a crowd, which seems to have sprung out of the earth, is gathered round the timber-yard.

    "It looks like a row, your honour..." says the policeman.

    Otchumyelov makes a half turn to the left and strides towards the crowd.

    He sees the aforementioned man in the unbuttoned waistcoat standing close by the gate of the timber-yard, holding his right hand in the air and displaying a bleeding finger to the crowd. On his half-drunken face there is a look of terror, and the finger itself has the look of a flag of victory. In this man Otchumyelov recognizes Hryukin, the goldsmith. The culprit who has caused the sensation, a white borzoy puppy with a sharp muzzle and a yellow patch on its back, is sitting on the ground with its front paws outstretched in the middle of the crowd, trembling all over. There is an expression of misery and terror in its tearful eyes.

    "What's it all about?" Otchumyelov inquires, pushing his way through the crowd. "What are you here for? Why are you holding up your finger?... Who was it shouted?"

    "I was walking along here, not interfering with anyone, your honour..." Hryukin begins, coughing into his fist. "I was talking about firewood with Mitry Mitritch, when this low brute, for no rhyme or reason, bit my finger.... You must excuse me, I am a working man.... Mine is fine work. I must have damages, for I shan't be able to use this finger for a week, may be.... It's not even the law, your honour, that one should put up with it from a beast.... If everyone is going to be bitten, life won't be worth living...."

    "H'm.... Very good," says Otchumyelov sternly, coughing and twitching his eyebrows. "Very good. Whose dog is it? I won't let this pass! I'll teach them to let their dogs run all over the place! It's time these gentry were looked after, if they won't obey the regulations! When he's fined, the blackguard, I'll teach him what it means to keep dogs and such stray cattle! I'll give him a lesson! Yeldrin," cries the superintendent, addressing the policeman, "find out who owns the dog and draw up a report! And the dog must be strangled. Without delay! It's sure to be mad.... Whose dog is it, I ask?"

    "I fancy it's General Zhigalov's," says someone in the crowd.

    "General Zhigalov's, h'm.... Take my coat off, Yeldrin," said Otchumyelov. "It's terribly hot! It must be a sign of rain.... There's one thing I can't make out, how it came to bite you?" Otchumyelov says, turning to Hryukin. "Surely it couldn't reach your finger. It's a little dog, and you are a great hulking fellow! You must have scratched your finger with a nail, and then the idea struck you to get damages for it. We know your sort! I know you, you devils!"

    "He put a cigarette in its mouth, for a joke, and it flew at him.... He is a nonsensical fellow, your honour."

    "That's a lie, Squinteye! You didn't see. If you didn't see, what's the use of your telling lies? The General's a sensible man and understands what's what. If his dog is mangy or has bitten anyone, he will have it shot and not break the law. And if everyone is to be bitten, life won't be worth living."

    "The general's cook is coming, we'll ask him.... Hi, Prokhor! Come here, my dear man! Look at this dog.... Is it one of yours?"

    "What an idea! We have never had one like that."

    "There's no need to waste time asking," says Otchumyelov. "It's a stray dog! There's no need to stand here talking about it.... If I say it's a stray dog, it is a stray dog.... It must be destroyed, that's all about it."

    "And yet it may be the General's," said the police-superintendent, thinking aloud. "It's not written on its face.... I saw one like it the other day in his yard."

    "It is the General's, that's certain!" said a voice in the crowd.

    "H'm, help me on with my overcoat, Yeldrin, my lad... the wind's getting up.... I am cold.... You take it to the General's, and inquire there. Say I found it and sent it. And tell them not to let it out into the street.... It may be a valuable dog, and if every swine goes sticking a cigarette in its mouth, it will soon be ruined. A dog is a delicate animal.... And you, you blockhead, put your hand down! It's no use your displaying your fool of a finger. It's your own fault!..."

    "The General's cook is coming out of that house; we'll ask him."

    "Yes, you do."

    "Prokhor, come here."

    "The dog is not ours," said the cook, "it belongs to the General's brother who came the other day. Our master does not care for hounds. His brother is fond of them."

    "Do you mean to say His Excellency's brother is here? Vladimir Ivanitch?" inquired Otchumyelov, and his whole face was lighted up with a smile of ecstasy. "Well, I never! and I did not know. Has he come on a visit?"

    "Yes."

    "Well, I never! He has been missing his brother," said Otchumyelov. "And so that's his dog? Delighted to hear it. Take it.... It is not a bad pup.... A lively creature.... Snapped at this fellow's finger. Ha-ha-ha! Come, why are you trembling? R-r-r.... The rogue is angry.... A nice little pup."

    Prokhor called the dog and went away with it.... The crowd laughed at Hryukin.

    "I'll make you smart yet!" Otchumyelov threatened him, and wrapping himself in his greatcoat, goes on his way through the market-place.
    """

    # Split the text into sentences for easier searching.
    sentences = re.split(r'(?<=[.!?])\s+', story_text)

    # Key phrases that indicate a symbolic action with the coat
    key_phrases = [
        "take my coat off",
        "help me on with my overcoat",
        "wrapping himself in his greatcoat"
    ]
    
    found_descriptions = []
    for sentence in sentences:
        for phrase in key_phrases:
            # Check if a key phrase is in the sentence, and ensure it's not already added
            if phrase in sentence.lower() and sentence not in found_descriptions:
                found_descriptions.append(sentence.strip())

    print("Found the following descriptions of Otchumyelov's coat that symbolize his shifting mentality:")
    print("-" * 80)
    
    equation_parts = []
    for i, desc in enumerate(found_descriptions, 1):
        print(f"{i}. {desc}")
        equation_parts.append("1")

    total_count = len(found_descriptions)
    equation_str = " + ".join(equation_parts)
    
    print("-" * 80)
    print(f"The total count is derived from each description: {equation_str} = {total_count}")

find_coat_descriptions()
>>>3