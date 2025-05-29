def check_adjective_order(sentence):
    adjectives = sentence.split()[:-1]  # Get all words except the last one (the noun)
    order = ["quantity", "quality", "size", "age", "shape", "color", "proper", "purpose"]
    adjective_types = {"old": "age", "drinking": "purpose"}
    
    indices = [order.index(adjective_types[adj]) for adj in adjectives if adj in adjective_types]
    return indices == sorted(indices)

sentence_A = "drinking old monkey"
sentence_B = "old drinking monkey"

result_A = check_adjective_order(sentence_A)
result_B = check_adjective_order(sentence_B)

print("Sentence A is correct:", result_A)
print("Sentence B is correct:", result_B)