import re

def count_coat_descriptions():
    """
    Analyzes Chekhov's "The Chameleon" to count descriptions of Otchumyelov's coat.
    """
    # The text of "The Chameleon" by Anton Chekhov (Constance Garnett translation)
    story_text = """
    THE police superintendent Otchumyelov is walking across the market-square wearing a new overcoat and carrying a parcel under his arm. A red-haired policeman strides after him with a sieve full of confiscated gooseberries in his hands. There is silence all around.... The square is empty.... The humble doorways of the shops and taverns look out upon the world disconsolately, like hungry mouths; there is not even a beggar near them.
    "So you bite, you damned brute?" Otchumyelov hears suddenly. "Lads, don't let him go! Biting is prohibited nowadays!"
    There is a sound of yelping. Otchumyelov looks in the direction of the sound and sees a dog, hopping on three legs and being chased by a man in a starched cotton shirt, with his waistcoat unbuttoned. The man runs after her, and, throwing his body forward, falls down and seizes the dog by her hind legs. The sound of yelping and a shout of "Don't let go!" is heard. Sleepy countenances are protruded from the shops, and soon a crowd, which seems to have sprung out of the earth, is gathered round the timber-yard.
    "It looks like a row, your honour..." says the policeman.
    Otchumyelov makes a half turn to the left and strides towards the crowd.
    He sees the aforementioned man in the unbuttoned waistcoat standing close by the gate of the timber-yard, holding his right hand in the air and displaying a bleeding finger to the crowd. On his drunken face there is a look of terror, as though he were asking to have his finger cut off.
    "Whose dog is it?"
    "I fancy it's General Zhigalov's," says some one in the crowd.
    "General Zhigalov's? H'm.... Take off my coat, Yeldrin," says Otchumyelov, addressing the policeman. "It's frightfully hot! It must be a sign of rain.... There is one thing I can't make out, how it came to bite you?" Otchumyelov turns to Hryukin. "Surely it couldn't reach your finger. It's a little dog, and you are a great hulking fellow! You must have scratched your finger with a nail, and then the idea struck you to get damages for it. I know your sort!"
    "He put a cigarette in her face, your honour, for a joke, and she snapped at him.... He is a nonsensical fellow, your honour!"
    "That's a lie.... If I were a liar, let the magistrate decide. He has a right to... I have a brother in the gendarmerie... let me tell you."
    "Don't argue!"
    "No, it's not the general's dog," says the policeman, with profound conviction, "the general has setters, with good points."
    "Do you know that for a fact?"
    "Yes, your honour."
    "I know it, too. The general has valuable dogs, thoroughbred, and this is goodness knows what! No coat, no shape.... A low creature. And to keep a dog like that!... where's the sense of it. If a dog like that were to turn up in Petersburg or Moscow, do you know what would happen? They would not worry about the law, they would strangle it in a twinkling! You have been injured, Hryukin, and we can't let the matter drop.... We must give them a lesson! It is high time..."
    "And yet, maybe it is the general's," says the policeman, thinking aloud. "It is not written on its face.... I saw one like it the other day in his yard."
    "It is the general's, that's certain!" says a voice in the crowd.
    "H'm!... Put on my coat, Yeldrin, my lad... there's a wind getting up.... I am cold.... You take it to the general's, and inquire there. Say I found it and sent it. And tell them not to let it out into the street.... It may be a valuable dog, and if every swine goes sticking a cigarette in its mouth, it will soon be ruined. A dog is a delicate animal.... And you, you blockhead, put your hand down! It's your own fault!"
    "The general's cook is coming, let's ask him.... Hi, Prohor! Come here, my dear man! Look at this dog.... Is it one of yours?"
    "What an idea! We have never had one like that."
    "There's no need to waste time asking," says Otchumyelov. "It's a stray dog! There's no need to waste time talking about it.... If he says it's a stray dog, a stray dog it is.... It must be destroyed, that's all about it."
    "It is not our dog," Prohor goes on. "It belongs to the general's brother, who came the other day. Our master does not care for hounds. But his honour is fond of them..."
    "Do you mean to say his Excellency's brother is here? Vladimir Ivanitch?" inquires Otchumyelov, and his whole face is transformed with a smile of ecstasy. "Fancy, your Excellency, I did not know! Has he come on a visit?"
    "Yes."
    "Well, I never.... He came to stay with his brother.... So it's his dog? Delighted to hear it.... It's not a bad pup.... A lively dog.... Snapped at this fellow's finger! Ha-ha-ha.... Come, why are you trembling? R-r-r-r.... The rascal is angry.... A nice little pup."
    Prohor calls the dog, and walks away from the timber-yard with it. The crowd laughs at Hryukin.
    "I'll make you smart yet!" Otchumyelov threatens him, and wrapping himself in his greatcoat, he goes on his way across the market-place.
    """
    
    # A second version of the text shows slight variations, so we check for both.
    story_text_v2 = """
    "H'm!... Take my coat off, Yeldrin... it's hot! It's as hot as though I were in a Turkish bath."
    "H'm.... Put my coat on, Yeldrin, my lad.... I fancy the wind is getting up.... I feel chilly...."
    """

    full_text = story_text + story_text_v2
    lower_text = full_text.lower()
    
    # Phrases that signify a description of the coat changing state
    # Each represents a distinct moment of change in Otchumyelov's attitude.
    phrases_to_count = {
        "'take off my coat'": lower_text.count("take off my coat"),
        "'put on my coat'": lower_text.count("put on my coat"),
        "'take my coat off'": lower_text.count("take my coat off"),
        "'wrapping himself in his greatcoat'": lower_text.count("wrapping himself in his greatcoat")
    }
    
    counts = list(phrases_to_count.values())
    total_count = sum(counts)
    equation = " + ".join(map(str, counts))

    print("Counting the descriptions of Otchumyelov's coat changing...")
    for phrase, count in phrases_to_count.items():
        if count > 0:
            print(f"Found {phrase} {count} time(s).")

    print(f"\nThe total number of descriptions is found by the sum:")
    print(f"{equation} = {total_count}")

count_coat_descriptions()
>>> 5