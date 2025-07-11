import collections

# Define a named tuple for better data organization
Painting = collections.namedtuple('Painting', ['label', 'title', 'artist', 'nationality_at_creation', 'year'])

# Store the information about each painting in a list
paintings_data = [
    Painting('A', 'Bathers on a Summer Evening', 'Félix Vallotton', 'Swiss', 1892),
    Painting('B', 'The Art Dealer (study)', 'Édouard Vuillard', 'French', 1898),
    Painting('C', 'The Olive Trees', 'Henri-Edmond Cross', 'French', 1892),
    Painting('D', 'Closed Eyes', 'Odilon Redon', 'French', 1890),
    Painting('E', 'The Hangover (Gueule de bois)', 'Henri de Toulouse-Lautrec', 'French', 1889),
    Painting('F', 'A Quiet Countryside', 'AI-generated/Modern Digital Art', 'N/A', 21) # 21st century
]

result_string = ""
print("Analyzing which paintings were created by a French painter before 1900:\n")

for p in paintings_data:
    # Condition 1: Is the artist French?
    is_french = p.nationality_at_creation == 'French'
    
    # Condition 2: Was the creation year before 1900?
    before_1900 = p.year < 1900

    # Check if both conditions are met
    if is_french and before_1900:
        decision = 'Y'
        reason = f"({p.artist} was a French painter and the work was created in {p.year}, which is before 1900.)"
    else:
        decision = 'N'
        reason = f"({p.artist} was {p.nationality_at_creation}, not French at the time of creation in {p.year}, or the work was created after 1900.)"
        if p.label == 'A':
            reason = f"({p.artist} was {p.nationality_at_creation} when the work was created in {p.year}. He became a French citizen in 1900.)"
        elif p.label == 'F':
            reason = "(This image is a modern digital creation, not a 19th-century painting.)"


    result_string += decision
    print(f"Painting {p.label}: {decision} {reason}")

print("\nFinal Answer String:")
print(result_string)

<<<NYYYYN>>>