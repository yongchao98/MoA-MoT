def find_true_sentences():
    """
    This function analyzes the components of four sentences to find all possible
    "true and meaningful" combinations by pairing subjects and objects within
    the original sentence structures.
    """
    subjects = ["Satyendranath Bose", "The conversation", "Sushi", "Our Tuna"]
    objects = ["Bosons", "playful", "fish", "a Bengalese"]
    
    true_sentences = []
    
    print("Finding true and meaningful sentences...\n")
    
    # Rule 1: Sentences with "baptized"
    # Semantic rule: A person (metaphorically) baptizes/names a concept.
    count1 = 0
    s = "Satyendranath Bose"
    o = "Bosons"
    sentence = f"{s} baptized {o}."
    true_sentences.append(sentence)
    count1 += 1
    print(sentence)
    
    # Rule 2: Sentences with "turned"
    # Semantic rule: An entity changes its state to a quality (adjective).
    # 'playful' is the only quality/adjective in the object list.
    count2 = 0
    o = "playful"
    for s in ["Satyendranath Bose", "The conversation", "Our Tuna"]:
        # "Sushi" is inanimate and cannot turn playful.
        sentence = f"{s} turned {o}."
        true_sentences.append(sentence)
        count2 += 1
        print(sentence)

    # Rule 3: Sentences with "contains"
    # Semantic rule: A composite item/being contains a physical component.
    # 'fish' is the component. 'Sushi' (a dish) and 'Our Tuna' (an animal that eats) can contain it.
    count3 = 0
    o = "fish"
    for s in ["Sushi", "Our Tuna"]:
        sentence = f"{s} contains {o}."
        true_sentences.append(sentence)
        count3 += 1
        print(sentence)
        
    # Rule 4: Sentences with "is a"
    # Semantic rule: A specific being is a member of a category.
    # "Our Tuna" (a cat) is a "Bengalese" (a cat breed). Other combinations are false.
    count4 = 0
    s = "Our Tuna"
    o = "a Bengalese"
    sentence = f"{s} is {o}."
    true_sentences.append(sentence)
    count4 += 1
    print(sentence)
    
    print("\n-------------------------------------------")
    total_count = len(true_sentences)
    print(f"The total number of true and meaningful sentences is the sum of the counts from each verb structure:")
    print(f"{count1} (from 'baptized') + {count2} (from 'turned') + {count3} (from 'contains') + {count4} (from 'is a') = {total_count}")

find_true_sentences()

<<<7>>>