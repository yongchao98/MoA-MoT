import collections

# A data structure to hold the analysis for each plot.
# Each entry is: Plot Number, (b, c, d) characteristics, parameters, code.
analysis = [
    (1, "Symmetric, Chaotic, Unbiased", "b=0, c=1, d=0", "C"),
    (2, "Asymmetric (red>blue), Medium, Unbiased", "b=1, c=0, d=0", "B"),
    (3, "Asymmetric (red>blue), Medium, Red bias", "b=1, c=0, d=1", "Z"),
    (4, "Asymmetric (blue>red), Chaotic, Blue bias", "b=-1, c=1, d=-1", "C"),
    (5, "Symmetric, Medium, Blue bias", "b=0, c=0, d=-1", "d"),
    (6, "Symmetric, Medium, Red bias", "b=0, c=0, d=1", "D"),
    (7, "Symmetric, Chaotic, Red bias", "b=0, c=1, d=1", "z"),
    (8, "Asymmetric (blue>red), Chaotic, Blue bias", "b=-1, c=1, d=-1", "C"),
    (9, "Asymmetric (blue>red), Medium, Unbiased", "b=-1, c=0, d=0", "b"),
    (10, "Symmetric, Ordered, Unbiased", "b=0, c=-1, d=0", "c"),
    (11, "Asymmetric (red>blue), Chaotic, Blue bias", "b=1, c=1, d=-1", "d"),
    (12, "Symmetric, Ordered, Blue bias", "b=0, c=-1, d=-1", "z"),
    (13, "Symmetric, Ordered, Blue bias", "b=0, c=-1, d=-1", "z"),
    (14, "Asymmetric (blue>red), Medium, Blue bias", "b=-1, c=0, d=-1", "Z"),
    (15, "Asymmetric (red>blue), Ordered, Red bias", "b=1, c=-1, d=1", "c"),
    (16, "Symmetric, Chaotic, Unbiased", "b=0, c=1, d=0", "C")
]

# Ensure the list is sorted by plot number, though it already is.
analysis.sort(key=lambda x: x[0])

# Extract the codes in order
result_string = "".join([item[3] for item in analysis])

print(result_string)
