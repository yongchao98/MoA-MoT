def check_adjective_order(sentence):
    order = ["quality", "age", "proper", "material"]
    adjectives = {
        "nice": "quality",
        "old-fashioned": "age",
        "Iranian": "proper",
        "lead": "material"
    }
    words = sentence.split()
    order_in_sentence = [adjectives[word] for word in words if word in adjectives]
    return order_in_sentence == order

sentence_A = "Iranian nice lead old-fashioned monkey"
sentence_B = "nice old-fashioned Iranian lead monkey"

result_A = check_adjective_order(sentence_A)
result_B = check_adjective_order(sentence_B)

print("Sentence A correct order:", result_A)
print("Sentence B correct order:", result_B)