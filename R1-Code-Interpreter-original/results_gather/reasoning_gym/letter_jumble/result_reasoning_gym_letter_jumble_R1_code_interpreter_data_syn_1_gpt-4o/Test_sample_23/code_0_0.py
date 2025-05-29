def unscramble(scrambled, possible):
    return sorted(scrambled) == sorted(possible)

words = [
    ("etxenivse", "extensive"),
    ("ni", "in"),
    ("eth", "the"),
    ("wdrol", "world"),
    ("cegishtrnt", "stretching")
]

results = [unscramble(scrambled, possible) for scrambled, possible in words]
print(results)