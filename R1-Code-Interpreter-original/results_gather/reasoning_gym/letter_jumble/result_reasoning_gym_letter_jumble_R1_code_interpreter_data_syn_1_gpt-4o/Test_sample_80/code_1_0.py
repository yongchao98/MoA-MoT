# Manually defined dictionary of possible words
word_dict = {
    "tse": "set",
    "orf": "for",
    "shi": "his",
    "reoirtncuesr": "reconstruction",
    "ti": "it",
    "swa": "was",
    "poposrde": "proposed",
    "ot": "to",
    "rM": "Mr",
    "hSmti": "Smith",
    "ttah": "that",
    "eh": "he",
    "ohsudl": "should",
    "rpmiet": "permit",
    "eth": "the",
    "snecdo": "second",
    "tpar": "part",
    "fo": "of"
}

# List of scrambled words
scrambled_words = ["tse", "orf", "shi", "reoirtncuesr", "ti", "swa", "poposrde", "ot", "rM", "hSmti", "ttah", "eh", "ohsudl", "rpmiet", "eth", "snecdo", "tpar", "fo"]

# Unscramble each word using the dictionary
unscrambled_sentence = ' '.join(word_dict.get(word, word) for word in scrambled_words)

print(unscrambled_sentence)