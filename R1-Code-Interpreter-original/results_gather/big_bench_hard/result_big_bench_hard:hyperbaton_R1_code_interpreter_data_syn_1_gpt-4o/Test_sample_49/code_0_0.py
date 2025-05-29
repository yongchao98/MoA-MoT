def check_adjective_order(sentence):
    order = ["size", "age", "shape", "proper adjective", "material"]
    adjectives = {
        "enormous": "size",
        "new": "age",
        "triangular": "shape",
        "Nigerian": "proper adjective",
        "lead": "material"
    }
    
    words = sentence.split()
    adjective_order = [adjectives[word] for word in words if word in adjectives]
    
    return adjective_order == order

sentence_A = "enormous new triangular Nigerian lead cat"
sentence_B = "enormous Nigerian triangular new lead cat"

result_A = check_adjective_order(sentence_A)
result_B = check_adjective_order(sentence_B)

print("Sentence A correct order:", result_A)
print("Sentence B correct order:", result_B)