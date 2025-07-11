def solve_riddle():
    """
    Solves the riddle by identifying a phonetic pun between Russian and English.
    """
    print("Step 1: The 'wide boulevards of Vienna' almost certainly refers to the famous Ringstrasse ('Ring Boulevard').")
    print("The key word here is 'Ring'.")
    print("-" * 20)

    print("Step 2: The riddle mentions a 'Russian translation'. We need the Russian word for 'Ring'.")
    russian_ring = "Кольцо"
    pronunciation_ring = "kol'tso"
    print(f"The Russian word for 'Ring' is '{russian_ring}', pronounced roughly as '{pronunciation_ring}'.")
    print("-" * 20)

    print("Step 3: The riddle is a pun. We are looking for an English poet whose name sounds like the Russian word for 'Ring'.")
    poet_name_english = "Coleridge"
    poet_name_russian = "Кольридж"
    pronunciation_poet = "Kol'ridzh"
    print(f"The surname of the English poet Samuel Taylor Coleridge is transliterated into Russian as '{poet_name_russian}' ({pronunciation_poet}).")
    print("-" * 20)

    print("Step 4: The first part of the surname, 'Коль-' (Kol'-), is a clear phonetic pun on 'Кольцо' (kol'tso).")
    print("Thus, in a Russian text about the Vienna 'Ring' (Кольцо), the name 'Coleridge' (Кольридж) is evoked.")
    print("-" * 20)
    
    final_answer = poet_name_english
    print(f"The surname is: {final_answer}")

solve_riddle()