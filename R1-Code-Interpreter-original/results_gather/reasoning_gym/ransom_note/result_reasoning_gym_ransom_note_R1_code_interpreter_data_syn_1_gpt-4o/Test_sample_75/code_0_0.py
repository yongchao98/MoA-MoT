from collections import Counter

def can_construct_ransom_note_from_magazine(ransom_note, magazine):
    ransom_note_count = Counter(ransom_note)
    magazine_count = Counter(magazine)
    
    for letter, count in ransom_note_count.items():
        if magazine_count[letter] < count:
            return False
    return True

ransom_note = "buhy"
magazine = "khsmjxvarzulyyxjcnqwyum"

result = can_construct_ransom_note_from_magazine(ransom_note, magazine)
print(result)